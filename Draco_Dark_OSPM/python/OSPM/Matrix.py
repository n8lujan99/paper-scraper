###################################################################################
# Matrix assembly and orbit weight solving routines for OSPM
###################################################################################

# Space grid setup
def spacel():
    pass

def sos2nuclei(Norbit, Nsos=55, nquadlength=10000, dist_lim=0.1, fracshift=0.001, minshift=1e-7, maxshift=1e5, tooclose=1e-8, do_plot_elz=False):
    def read_sos():
        """Read sos1.out and sos2.out -> rs, vs, ts arrays + bins."""
        mlnorm = float(np.loadtxt("gden.norm"))
        df2 = pd.read_csv("sos2.out", delim_whitespace=True, names=["iorb", "iE", "iLz", "dLz", "dE"])
        df1 = pd.read_csv("sos1.out", delim_whitespace=True, names=["rs_in", "vs_in", "ts_in"])

        fac1 = 1.0  # external scaling not applied here
        rmax = 1.0

        rs = np.full((Norbit, Nsos), np.nan)
        vs = np.full_like(rs, np.nan)
        ts = np.full_like(rs, np.nan)
        E_bin = np.zeros(Norbit, int)
        Lz_bin = np.zeros(Norbit, int)
        nsosorb = np.zeros(Norbit, int)

        for i, (iorb, iE, iLz, dLz, dE) in enumerate(df2.values):
            rs_in, vs_in, ts_in = df1.iloc[i]
            iorb = int(iorb) - 1  # Fortran 1-based -> Python 0-based
            nsosorb[iorb] += 1
            idx = nsosorb[iorb] - 1
            rs[iorb, idx] = np.log10(rs_in * rmax)
            vs[iorb, idx] = vs_in * fac1 * mlnorm
            ts[iorb, idx] = ts_in
            E_bin[iorb] = int(iE)
            Lz_bin[iorb] = int(iLz)

        return rs, vs, ts, nsosorb, E_bin, Lz_bin

    # === helper: group orbits into unique (E,Lz) sequences =====================
    def sosseq(E_bin, Lz_bin):
        data = pd.DataFrame({"E": E_bin, "Lz": Lz_bin})
        grouped = data.drop_duplicates().reset_index(drop=True)
        E_seq = grouped["E"].values
        Lz_seq = grouped["Lz"].values
        nseq = len(E_seq)
        seqoforb = np.zeros(len(E_bin), int)
        for iorb, (E, Lz) in enumerate(zip(E_bin, Lz_bin)):
            seqoforb[iorb] = np.where((E_seq == E) & (Lz_seq == Lz))[0][0]
        return nseq, E_seq, Lz_seq, seqoforb

    # === helper: plot ELz like Fortran plotELz (optional) =====================
    def _plot_ELz():
        """Replicates plotELz: scatter of log10(-E) vs ±Lz from integrals.out."""
        # integrals.out lines: xLz, xE, x3, x4; Fortran used E=log10(-x2), Lz=±x1
        arr = np.loadtxt("integrals.out")
        xLz = arr[:, 0]
        xE  = arr[:, 1]
        E_vals = np.log10(-xE)
        Lz_vals = xLz
        # duplicate for symmetric ±Lz, as in Fortran
        E_plot = np.repeat(E_vals, 2)
        Lz_plot = np.empty(2 * len(Lz_vals))
        Lz_plot[0::2] = -Lz_vals
        Lz_plot[1::2] =  Lz_vals

        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(E_plot, Lz_plot, s=10)
        ax.set_xlabel(r'$\log_{10}(-E)$')
        ax.set_ylabel(r'$L_z$')
        ax.set_title('ELz diagram (integrals.out)')
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        return fig, ax

    # ===================== main routine =======================================
    rs, vs, ts, nsosorb, E_bin, Lz_bin = read_sos()
    nseq, E_seq, Lz_seq, seqoforb = sosseq(E_bin, Lz_bin)

    nuclei_list, vor_cells_list, vor_sos_list, vor_plot_list = [], [], [], []
    seq_summary = []

    for iseq in range(nseq):
        # collect all SOS points for this sequence
        idxs = np.where(seqoforb == iseq)[0]
        x, y, tarr, orbarr = [], [], [], []
        for iorb in idxs:
            count = int(nsosorb[iorb])
            if count <= 0:
                continue
            x.extend(10 ** rs[iorb, :count])
            y.extend(vs[iorb, :count])
            tarr.extend(ts[iorb, :count])
            orbarr.extend([iorb + 1] * count)

        if not x:
            continue

        x = np.asarray(x); y = np.asarray(y)
        tarr = np.asarray(tarr); orbarr = np.asarray(orbarr)

        # normalize x,y by their maxima
        xmax, ymax = np.max(x), np.max(y)
        x = x / xmax
        y = y / ymax

        # --- build envelope (upper ridge), then mirror across v=0
        i_max = np.argmax(y)
        xt = [x[i_max]]
        yt = [y[i_max] + np.clip(fracshift * y[i_max], minshift, maxshift)]

        # to larger radii
        xcur = x[i_max]
        while True:
            idx = np.where(x > xcur)[0]
            if not len(idx):
                break
            i_next = idx[np.argmax(y[idx])]
            if i_next == i_max:
                break
            xcur = x[i_next]
            xt.append(xcur + np.clip(fracshift * xcur, minshift, maxshift))
            yt.append(y[i_next] + np.clip(fracshift * y[i_next], minshift, maxshift))
            if np.isclose(xcur, x.max()):
                break

        # to smaller radii
        xcur = x[i_max]
        while True:
            idx = np.where(x < xcur)[0]
            if not len(idx):
                break
            i_next = idx[np.argmax(y[idx])]
            if i_next == i_max:
                break
            xcur = x[i_next]
            xt.append(xcur * (1 - fracshift))
            yt.append(y[i_next] * (1 + fracshift))
            if np.isclose(xcur, x.min()):
                break

        xt = np.asarray(xt); yt = np.asarray(yt)
        order = np.argsort(xt)
        xt, yt = xt[order], yt[order]

        # interpolate to densify boundary (Fortran's ysmall logic)
        ysmall = max(min((xt[-1] - xt[0]) / 40.0, 1e-2), 1e-3)
        xti, yti = [xt[0]], [yt[0]]
        for i in range(1, len(xt)):
            dx = xt[i] - xt[i - 1]; dy = yt[i] - yt[i - 1]
            steps = int(max(abs(dx), abs(dy)) / ysmall)
            for j in range(1, steps):
                f = j / steps
                xti.append(xt[i - 1] + f * dx)
                yti.append(yt[i - 1] + f * dy)
            xti.append(xt[i]); yti.append(yt[i])
        xt = np.asarray(xti); yt = np.asarray(yti)

        # mirror points with y<dist_lim
        xenv = np.concatenate([1.001 * x[y < dist_lim], xt])
        yenv = np.concatenate([-y[y < dist_lim],       yt])

        # back to physical scale
        x_all = xenv * xmax
        y_all = yenv * ymax
        orb_all = np.concatenate([np.zeros_like(xenv, int), np.ones_like(xt, int)])
        t_all = np.zeros_like(x_all)

        # drop pairs that are too close (scaled by nquadlength like Fortran)
        valid = np.ones(len(x_all), dtype=bool)
        for i in range(len(x_all)):
            if not valid[i]:
                continue
            dx = np.abs(x_all[i] - x_all[i+1:]) / float(nquadlength)
            dy = np.abs(y_all[i] - y_all[i+1:]) / float(nquadlength)
            bad = (dx < tooclose) & (dy < tooclose)
            valid[i+1:][bad] = False

        x_all, y_all, orb_all, t_all = x_all[valid], y_all[valid], orb_all[valid], t_all[valid]

        # assemble outputs analogous to files
        N = len(x_all)
        xjunk = np.full(N, 0.30)
        indices = np.arange(N)
        nuclei   = np.column_stack((x_all, y_all, xjunk, indices))
        vor_cells= np.tile([1, 1, 1, -1.0], (N, 1))
        vor_sos  = np.column_stack((x_all, y_all, xjunk, orb_all, t_all, np.full(N, xmax), np.full(N, ymax)))
        vor_plot = np.column_stack((x_all, y_all, orb_all, np.zeros(N)))

        nuclei_list.append(nuclei)
        vor_cells_list.append(vor_cells)
        vor_sos_list.append(vor_sos)
        vor_plot_list.append(vor_plot)
        seq_summary.append((iseq + 1, E_seq[iseq], Lz_seq[iseq], N))

    seq_df = pd.DataFrame(seq_summary, columns=["seq", "E_seq", "Lz_seq", "N_points"])

    out = {
        "nseq": nseq,
        "nuclei": nuclei_list,
        "vor_cells": vor_cells_list,
        "vor_sos": vor_sos_list,
        "vor_plot": vor_plot_list,
        "seq_info": seq_df,
    }

    if do_plot_elz:
        fig, ax = _plot_ELz()
        out["elz_fig"] = fig
        out["elz_ax"]  = ax

    return out

def velmtod(veld, velm, ad, adfer, xmu, Nvelb):
    """Finds which velocity bins of data correspond to model velocity bin edges."""
    nad = veld.shape[0]
    Nvel = velm.shape[0]

    sumad = np.zeros((Nvel, Nvelb))
    sadfer = np.zeros((Nvel, Nvelb))

    vbin = velm[1] - velm[0]

    # Loop over apertures
    for irc in range(Nvelb):
        # --- Compute data bin edges
        v1 = np.zeros(nad)
        v2 = np.zeros(nad)

        v1[0] = veld[0, irc] * xmu - (veld[1, irc] - veld[0, irc]) * xmu / 2.0
        v2[0] = veld[0, irc] * xmu + (veld[1, irc] - veld[0, irc]) * xmu / 2.0
        v1[-1] = veld[-1, irc] * xmu - (veld[-1, irc] - veld[-2, irc]) * xmu / 2.0
        v2[-1] = veld[-1, irc] * xmu + (veld[-1, irc] - veld[-2, irc]) * xmu / 2.0

        for k in range(1, nad - 1):
            v1[k] = veld[k, irc] * xmu - (veld[k, irc] - veld[k - 1, irc]) * xmu / 2.0
            v2[k] = veld[k, irc] * xmu + (veld[k + 1, irc] - veld[k, irc]) * xmu / 2.0

        sum1 = np.sum(ad[:, irc])  # total amplitude

        # --- Map data bins into model bins
        for i in range(Nvel):
            vpos1 = velm[i] - vbin / 2.0
            vpos2 = velm[i] + vbin / 2.0
            total = 0.0
            totale = 0.0

            for k in range(nad):
                frac = 0.0
                if vpos1 < v1[k] and vpos2 > v2[k]:
                    frac = 1.0
                elif vpos1 > v1[k] and vpos1 < v2[k]:
                    frac = (v2[k] - vpos1) / (v2[k] - v1[k])
                elif vpos2 > v1[k] and vpos2 < v2[k]:
                    frac = (vpos2 - v1[k]) / (v2[k] - v1[k])

                total += frac * ad[k, irc]
                totale += frac * adfer[k, irc]

            if velm[i] < v1[0] or velm[i] > v2[-1]:
                total = 0.0
                totale = -666.0

            sumad[i, irc] = total
            sadfer[i, irc] = totale

    return sumad, sadfer

###################################################################################
# Assemble Matrices
###################################################################################

def cmatrix (xlib, v1lib, Norbit, Norbitm, Nvel, Nvelb, Nvelbm, Nbin):
    # Assembles the matrix in Cm for the matrix Am in spear
    """Assembles the matrix in Cm for the matrix Am in spear."""
    lds = Norbitm + Nvel * Nvelbm
    lda = Nbin + Nvel * Nvelbm
    Cm = np.zeros((lda, lds), dtype=float)
    Nvtot = Nvel * Nvelb

    # Left two blocks (upper-left and lower-left)
    for i in range(Norbit):
        # upper-left piece
        for j in range(Nbin):
            Cm[j, i] = float(xlib[j, i])

        # lower-left piece
        for irc in range(Nvelb):
            for ivel in range(Nvel):
                j = Nbin + irc * Nvel + ivel
                Cm[j, i] = float(v1lib[ivel, irc, i])

    # Right two pieces (initialize zeros)
    for i in range(Norbit, Norbit + Nvtot):
        for j in range(Nbin + Nvel * Nvelb):
            Cm[j, i] = 0.0

    # Identity entries in lower-right block
    for i in range(Nvtot):
        Cm[Nbin + i, Norbit + i] = 1.0

    return Cm

def getlowright(Cm, w, sumad, Nvel, ircmax, Norbit, Norbitm, Nvelbm, Nbin): 
    """assembles the lower left portion of the matrix Cm which is used to assemble the matrix Am in SPEAR"""
    lds = Norbitm + Nvel * Nvelbm
    lda = Nbin + Nvel * Nvelbm
    big = 1.0e20

    Nvtot = Nvel * ircmax
    ilast = Norbit + Nvtot      # Fortran's +1 skipped since Python is 0-indexed
    xmu = w[ilast]

    for irc in range(ircmax):
        for ivel in range(Nvel):
            i = Nbin + irc * Nvel + ivel
            Cm[i, ilast] = float(-sumad[ivel, irc])

    return Cm 

def newcmatrix(Cm, xlib, v1lib, inew, Norbit, Nvel, Nvelb, Nbin, Nvelbm):
    """Updates matrix Cm by swapping in orbit 'inew' and reassembling blocks."""
    Nvtot = Nvel * Nvelb
    ilast = Norbit + Nvtot

    # Temporary storage
    xlibtemp = np.zeros(Nbin, dtype=float)
    v1libtemp = np.zeros((Nvel, Nvelb), dtype=float)

    # Copy the orbit column to temporary arrays
    for i in range(Nbin):
        xlibtemp[i] = xlib[i, inew]

    for i in range(Nvelb):
        for ivel in range(Nvel):
            v1libtemp[ivel, i] = v1lib[ivel, i, inew]

    # Swap selected orbit with the last one (Norbit)
    for i in range(Nbin):
        xlib[i, inew] = xlib[i, Norbit]
        xlib[i, Norbit] = xlibtemp[i]

    for i in range(Nvelb):
        for ivel in range(Nvel):
            v1lib[ivel, i, inew] = v1lib[ivel, i, Norbit]
            v1lib[ivel, i, Norbit] = v1libtemp[ivel, i]

    # Upper-left block (density)
    for j in range(Nbin):
        Cm[j, inew] = xlib[j, inew]

    # Lower-left block (velocity)
    for irc in range(Nvelb):
        for ivel in range(Nvel):
            j = Nbin + irc * Nvel + ivel
            Cm[j, inew] = v1lib[ivel, irc, inew]

    # Decrement Norbit (orbit count reduced by one)
    Norbit -= 1

    # Right-side blocks (zeroed region)
    for i in range(Norbit, Norbit + Nvtot):
        for j in range(Nbin + Nvel * Nvelb):
            Cm[j, i] = 0.0

    # Lower-right identity entries
    for i in range(Nvtot):
        Cm[Nbin + i, Norbit + i] = 1.0

    return Cm, xlib, v1lib, Norbit

def spear(Cm, S, dS, ddS, rcond, negweights, apfac2, Norbit, Nvel, Nvelb, Nbin, Nvelbm, w, xlib, v1lib, SMc, sumad, fracnew, apfac):
    """Compute Lagrange multipliers and Newton–Raphson corrections for orbit weights."""

    # --- Define matrix dimensions (same as Fortran parameter() section)
    lds = Norbit + Nvel * Nvelbm
    lda = Nbin + Nvel * Nvelbm
    Nvtot = Nvel * Nvelb
    Narr = Nbin + Nvtot      # total number of observable constraints
    ilast = Norbit + Nvtot   # total number of active weight elements

    # 1. Construct scaled version of Cm
    Cma = Cm[:Narr, :ilast] / ddS[:ilast]

    # 2. Compute curvature matrix Am = Cm * Cma^T
    Am = Cm[:Narr, :ilast] @ Cma.T

    # 3. Compute residual vector delY
    #     delY contains differences between observed data and model prediction
    #     First part (0:Nbin): surface density residuals
    delY = np.zeros(Narr)
    for ibin in range(Nbin):
        # sum over all orbits contributing to this spatial bin
        sum_val = np.sum(w[:Norbit] * xlib[ibin, :Norbit])
        delY[ibin] = SMc[ibin] - sum_val

    # 4. Second part (Nbin:): velocity distribution residuals
    for irc in range(Nvelb):
        for ivel in range(Nvel):
            i = Nbin + irc * Nvel + ivel
            # predicted velocity contribution from all orbits
            sum_val = np.sum(w[:Norbit] * v1lib[ivel, irc, :Norbit])
            # plus constant velocity term weights
            sum_val += w[Norbit + irc * Nvel + ivel]
            # observed - model difference, scaled by fracnew
            delY[i] = sumad[ivel, irc] * fracnew - sum_val

    # 5. Add the C*dS/ddS correction term
    delY += Cm[:Narr, :ilast] @ (dS[:ilast] / ddS[:ilast])

    # Solve Am * delY = delY  (Fortran's dsico + dsisl)
    # This finds the Lagrange multipliers λ
    lu, piv = lu_factor(Am)
    delY = lu_solve((lu, piv), delY)

    # Adjust orbit weights (Newton–Raphson correction)

    for j in range(ilast):
        sum_val = np.dot(delY, Cm[:Narr, j])
        dw = (sum_val - dS[j]) / ddS[j]
        wnew = w[j] + apfac * dw
        w[j] = wnew

    return w





