      PROGRAM library
      INCLUDE 'libdefs.h'
      real raban,telephon

C---- open a file to record function and subroutine calls:

      open (unit=50,file='progressL.out',status='unknown')
      write (50,*)''


C---- read in input files:

      write (50,*)' -> reading in data and parameter files'
      call binset(Nrdat,Nvdat,Nrlib,Nvlib)
      call dataread()
      call galaxyread()
      write (50,*)'   - data acquisition complete'
      write (50,*)''

C---- read DM-halo parameters
      call haloread()

C---- determine if density distribution is a sphere and set volume arrays

      write (50,*)' -> setting volume arrays and quadrupole flag'
      call spacel()
      write (50,*)'   - done'
      write (50,*)''


C---- calculate core mass (sans hole) and central density

      write (50,*)' -> calculating core quantities'
      cden = dL(1,1)*ratML(1,1) / vol2d(1)
      core = 4.*pi*cden*rmin**3/(3.+cslope)
      write (50,*)'   - done'
      write (50,*)''
      
C---- calculate density profile, total mass, total light; write to file:

      write (50,*)' -> calculating density distribution'
      call density()
      call densitywrite()
c      call halodensity()
      call totmlwrite()
      write (50,*)'   - mass and light computed and written to file'
      write (50,*)''


C---- calculate the force tables Tabth and Tabr:

C---- in case the depro is required outside the field of view of the library
      if(iextra.eq.-1) then
         call rhobeyond()
         write (50,*)'   - depro used outside fov of the library'
         write (50,*)''
      end if
         
C---- now proceed as normal
      write (50,*)' -> calculating force tables:'
C   pick one of either tablesread() or tables():
c      call tablesread()
      call tables()
      call halotables()
c      call halotableswrite()
      call addtables()
      call tableswrite()
      write (50,*)'   - force tables computed'
      write (50,*)''

C---- take care of non-spherical problems with core and outer edge

      if(iquad.ne.0) call forcefix()


      call potwrite()

C---- sample available orbits:

      write (50,*)' -> calling subroutine librarian'
      open (unit=79,file='integrals.out',status='unknown')
      open (unit=80,file='librarian.out',status='unknown')
      call librarian()
      close(79)
      close(80)
      write (50,*)'   - phase space sampled:  orbit library created'
      write (50,*)''


C---- write out Output files:

      write (50,*)' -> writing out phase space densities'
      call phasewrite()
      write (50,*)'   - phase space densities written'
      write (50,*)''

      open(unit=1,file='norbit',status='unknown')
      write(1,*) Norbtot
      close(1)

      close(50)

      END
