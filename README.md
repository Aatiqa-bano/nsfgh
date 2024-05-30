# nsfgh
import tkinter as tk
from tkinter import messagebox, filedialog

def validate_sequence(sequence, valid_chars):
    sequence = sequence.upper()  # Convert sequence to uppercase
    return all(base in valid_chars for base in sequence)

def transcribe_dna_to_rna(dna_sequence):
    return dna_sequence.replace('T', 'U')

def translate_sequence(sequence, is_rna):
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '', 'UAG': '',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    amino_acids = ""
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]

        # Check for stop codons
        if is_rna:
            stop_codons = ['UAA', 'UAG', 'UGA']
        else:
            stop_codons = ['TAA', 'TAG', 'TGA']

        if codon in stop_codons:
            break

        amino_acid = codon_table.get(codon, '')
        if amino_acid:
            amino_acids += amino_acid

    return amino_acids or None

def get_reading_frames(sequence, is_rna):
    frames = []
    start_index = 0

    while start_index < len(sequence):
        # Find the start codon 'ATG' or 'AUG' for RNA
        start_codon = 'AUG' if is_rna else 'ATG'
        start_index = sequence.find(start_codon, start_index)

        if start_index == -1:
            break  # If 'ATG' or 'AUG' is not found, exit the loop

        # Forward reading frame
        forward_frame = sequence[start_index:]
        frames.append(translate_sequence(forward_frame, is_rna))
        start_index += 3  # Move to the next position to avoid overlapping frames

    return frames[:3], frames[3:6]  # Return only the first 3 forward and reverse frames

def translate():
    sequence = entry.get()  # Do not convert input sequence to uppercase
    print("Input Sequence:", sequence)  # Debugging statement
    choice = var.get()

    valid_dna_chars = {'A', 'T', 'C', 'G'}
    valid_rna_chars = {'A', 'U', 'C', 'G'}

    print("Choice:", choice)  # Debugging statement

    if sequence:  # Check if sequence is provided directly
        if choice == 'DNA':
            if validate_sequence(sequence, valid_dna_chars):
                print("Valid DNA sequence.")  # Debugging statement
                forward_frames, reverse_frames = get_reading_frames(sequence, False)
                display_output(forward_frames, reverse_frames)
            else:
                messagebox.showerror("Error", "Invalid DNA sequence. Please enter a valid DNA sequence.")
        elif choice == 'RNA':
            if validate_sequence(sequence, valid_rna_chars):
                print("Valid RNA sequence.")  # Debugging statement
                forward_frames, reverse_frames = get_reading_frames(sequence, True)
                display_output(forward_frames, reverse_frames)
            else:
                messagebox.showerror("Error", "Invalid RNA sequence. Please enter a valid RNA sequence.")
        else:
            messagebox.showerror("Error", "Invalid choice. Please select DNA or RNA.")
    else:
        messagebox.showerror("Error", "Please enter a sequence.")

def display_output(forward_frames, reverse_frames):
    # Clear previous output
    output_text.delete(1.0, tk.END)

    # Print amino acid sequences for each reading frame
    for i, frame in enumerate(forward_frames):
        if frame:
            output_text.insert(tk.END, f"Amino acid sequence (Forward Frame {i + 1}): {frame.replace('_', '')}\n")
        else:
            output_text.insert(tk.END, f"No amino acid sequence found for Forward Frame {i + 1}\n")

    for i, frame in enumerate(reverse_frames):
        if frame:
            output_text.insert(tk.END, f"Amino acid sequence (Reverse Frame {i + 1}): {frame.replace('_', '')}\n")
        else:
            output_text.insert(tk.END, f"No amino acid sequence found for Reverse Frame {i + 1}\n")

def reset():
    entry.delete(0, tk.END)
    output_text.delete(1.0, tk.END)

def browse_file_translation():
    file_path = filedialog.askopenfilename(filetypes=[("FASTA files", ".fasta"), ("All files", ".*")])
    if file_path:
        with open(file_path, 'r') as file:
            # Read the contents of the file
            lines = file.readlines()

            # Concatenate sequence lines, removing header lines and white spaces
            sequence_lines = [line.strip() for line in lines if not line.startswith('>') and not line.isspace()]
            sequence = ''.join(sequence_lines)

            # Add the sequence to the entry field
            entry.delete(0, tk.END)
            entry.insert(0, sequence)     


def save_results():
    # Check if the current operation is forward translation
    choice = var.get()  # Get the choice for forward translation

    if choice in ('DNA', 'RNA'):
        # Forward translation
        result_widget = output_text
        file_extension = '.txt'  # Default file extension for forward translation
        file_types = [("Text files", ".txt"), ("All files", ".*")]

    # Check if there's any result to save
    result = result_widget.get(1.0, tk.END)
    if not result.strip():
        messagebox.showwarning("No Result", "No translation result to save.")
        return

    # Prompt user to choose the save location and file name
    file_path = filedialog.asksaveasfilename(defaultextension=file_extension, filetypes=file_types)
    if file_path:
        try:
            # Write the result to the selected file
            with open(file_path, 'w') as file:
                file.write(result)
            messagebox.showinfo("Saved", "Translation result saved successfully.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while saving the translation result: {e}")
            print(f"Error occurred: {e}")

# Create window
window = tk.Tk()
window.title("DNA/RNA Translator")

# Create label for input
input_label = tk.Label(window, text="Enter DNA/RNA sequence or load from file:", font=("Arial", 12))
input_label.grid(row=0, column=0, padx=10, pady=5, sticky='w')

# Create entry for sequence input
entry = tk.Entry(window, width=50, font=("Arial", 12))
entry.grid(row=1, column=0, padx=10, pady=5, sticky='w')

# Create browse button for translation
browse_button_translation = tk.Button(window, text="Browse the Fasta file", command=browse_file_translation, font=("Arial", 12))
browse_button_translation.grid(row=2, column=0, padx=10, pady=5, sticky='w')

# Create radio buttons for selecting DNA or RNA
var = tk.StringVar(value="DNA")
dna_radio = tk.Radiobutton(window, text="DNA", variable=var, value="DNA", font=("Arial", 10))
dna_radio.grid(row=1, column=1, padx=10, pady=5, sticky='w')
rna_radio = tk.Radiobutton(window, text="RNA", variable=var, value="RNA", font=("Arial", 10))
rna_radio.grid(row=1, column=2, padx=10, pady=5, sticky='w')

# Create translate button
translate_button = tk.Button(window, text="Translate", command=translate, font=("Arial", 12), bg="#4CAF50", fg="white", relief=tk.FLAT, padx=10)
translate_button.grid(row=3, column=0, padx=10, pady=5, sticky='w')

# Bind the Enter key to the translate function
window.bind('<Return>', lambda event=None: translate())

# Create reset button
reset_button = tk.Button(window, text="Reset", command=reset, font=("Arial", 12), bg="#f44336", fg="white", relief=tk.FLAT, padx=10)
reset_button.grid(row=3, column=1, padx=10, pady=5, sticky='w')

# Create save button
save_button = tk.Button(window, text="Save Results", command=save_results, font=("Arial", 12))
save_button.grid(row=3, column=2, padx=10, pady=5, sticky='w')

# Create text widget for output
output_text = tk.Text(window, height=10, width=70, font=("Arial", 12))
output_text.grid(row=4, column=0, padx=10, pady=10, columnspan=3)

# Run the GUI
window.mainloop()
