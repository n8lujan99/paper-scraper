import itertools
import numpy as np
import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog
import os


def generate_combinations(ml_range, bh_mass_range, dm_vcirc_range, dm_rscale_range, kinematics_value):
    for comb in itertools.product(ml_range, bh_mass_range, dm_vcirc_range, dm_rscale_range):
        comb = tuple("{:.0f}".format(value) if i == 1 else value for i, value in enumerate(comb))
        combination_str = ' '.join(format(value, ".1f") if isinstance(value, float) else str(value) for value in comb) + f' {kinematics_value}'
        yield "runallt " + combination_str


def write_combinations(start_index, end_index, file_index, ml_range, bh_mass_range, dm_vcirc_range, dm_rscale_range, kinematics_value, save_dir):
    filename = os.path.join(save_dir, f"combinations_part_{file_index}.txt")
    with open(filename, "w") as file:
        generator = generate_combinations(ml_range, bh_mass_range, dm_vcirc_range, dm_rscale_range, kinematics_value)
        for _ in range(start_index):
            next(generator)  # Skip to the start index
        for _ in range(start_index, end_index):
            file.write(next(generator) + '\n')
    print(f"File created: {filename}")


def on_submit():
    try:
        ml_start = float(entry_ml_start.get())
        ml_end = float(entry_ml_end.get())
        ml_step = float(entry_ml_step.get())

        bh_mass_start = float(entry_bh_mass_start.get())
        bh_mass_end = float(entry_bh_mass_end.get())
        bh_mass_step = float(entry_bh_mass_step.get())

        dm_vcirc_start = float(entry_dm_vcirc_start.get())
        dm_vcirc_end = float(entry_dm_vcirc_end.get())
        dm_vcirc_step = float(entry_dm_vcirc_step.get())

        dm_rscale_start = float(entry_dm_rscale_start.get())
        dm_rscale_end = float(entry_dm_rscale_end.get())
        dm_rscale_step = float(entry_dm_rscale_step.get())

        kinematics_value = int(entry_kinematics_value.get())
        if kinematics_value not in [0, 1]:
            raise ValueError("Kinematics value must be 0 or 1")

        ml_range = np.arange(ml_start, ml_end + ml_step, ml_step)
        bh_mass_range = np.arange(bh_mass_start, bh_mass_end + bh_mass_step, bh_mass_step)
        dm_vcirc_range = np.arange(dm_vcirc_start, dm_vcirc_end + dm_vcirc_step, dm_vcirc_step)
        dm_rscale_range = np.arange(dm_rscale_start, dm_rscale_end + dm_rscale_step, dm_rscale_step)

        total_combinations = len(ml_range) * len(bh_mass_range) * len(dm_vcirc_range) * len(dm_rscale_range)

        messagebox.showinfo("Result", f"Total Combinations: {total_combinations}")

        if messagebox.askyesno("Generate Files", "Do you want to generate the combination files?"):
            prompt_message = f"Enter the number of files to divide the {total_combinations} combinations into:"
            num_files = simpledialog.askinteger("Number of Files", prompt_message)
            if num_files is None or num_files <= 0:
                messagebox.showerror("Input Error", "Please enter a valid number of files.")
                return

            combinations_per_part = total_combinations // num_files  # Divide into the specified number of parts

            save_dir = filedialog.askdirectory(title="Select Directory to Save Files")
            if save_dir:
                split_indices = [(i * combinations_per_part) for i in range(num_files)]  # Create indices for splitting
                split_indices.append(total_combinations)  # Add the last index

                for i in range(num_files):
                    write_combinations(split_indices[i], split_indices[i + 1], i + 1, ml_range, bh_mass_range, dm_vcirc_range, dm_rscale_range, kinematics_value, save_dir)
                messagebox.showinfo("Done", "Files created successfully.")
    except ValueError as e:
        messagebox.showerror("Input Error", str(e))


# Create the main window
root = tk.Tk()
root.title("Parameter Combination Generator")

# Create and place the input fields and labels
labels = ["ML Start:", "ML End:", "ML Step:",
          "BH Mass Start:", "BH Mass End:", "BH Mass Step:",
          "DM VCirc Start:", "DM VCirc End:", "DM VCirc Step:",
          "DM RScale Start:", "DM RScale End:", "DM RScale Step:",
          "Kinematics Value (0 or 1):"]
variables = [tk.StringVar() for _ in labels]
entries = []

for i, label_text in enumerate(labels):
    label = tk.Label(root, text=label_text)
    label.grid(row=i, column=0, sticky="e")
    entry = tk.Entry(root, textvariable=variables[i])
    entry.grid(row=i, column=1)
    entries.append(entry)

(entry_ml_start, entry_ml_end, entry_ml_step,
 entry_bh_mass_start, entry_bh_mass_end, entry_bh_mass_step,
 entry_dm_vcirc_start, entry_dm_vcirc_end, entry_dm_vcirc_step,
 entry_dm_rscale_start, entry_dm_rscale_end, entry_dm_rscale_step,
 entry_kinematics_value) = entries

# Set default values (optional)
defaults = ["0.2", "1.2", "0.1",
            "360000", "400000", "10000",
            "14.1", "14.9", "0.1",
            "0.5", "1", "0.1",
            "0"]  # Default kinematics value is set to 0
for var, default in zip(variables, defaults):
    var.set(default)

# Create and place the submit button
submit_button = tk.Button(root, text="Submit", command=on_submit)
submit_button.grid(row=len(labels), columnspan=2)

# Print the current working directory
print(f"Current working directory: {os.getcwd()}")

# Run the application
root.mainloop()
