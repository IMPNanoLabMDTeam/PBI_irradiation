import os
import sys

def clean_data_file(input_file, output_file):
    """
    Remove force field and topology information from LAMMPS data files.
    Keeps only header, atom masses, and atom coordinates.
    """
    # Read the file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Output lines
    output_lines = []
    
    # First line (title)
    output_lines.append(lines[0])
    output_lines.append("# This file has been cleaned by removing force field and topology information\n")
    
    # Process the header section
    i = 1
    while i < len(lines):
        line = lines[i].strip()
        # Add the number of atoms and atom types
        if "atoms" in line and not any(x in line for x in ["bonds", "angles", "dihedrals", "impropers"]):
            output_lines.append(lines[i])
        elif "atom types" in line:
            output_lines.append(lines[i])
        # Skip all bond types, angle types etc.
        elif any(x in line for x in ["bond types", "angle types", "dihedral types", "improper types"]):
            i += 1
            continue
        # Add the box dimensions
        elif any(x in line for x in ["xlo", "ylo", "zlo"]):
            output_lines.append(lines[i])
        # Stop when we reach the first section
        elif line.startswith("Masses") or line.startswith("Pair Coeffs"):
            break
        i += 1
    
    # Add an empty line before Masses
    output_lines.append("\n")
    
    # Find and add the Masses section
    masses_section = []
    in_masses = False
    i = 0
    
    while i < len(lines):
        if lines[i].strip() == "Masses":
            in_masses = True
            masses_section.append(lines[i])
            i += 1
            
            # Skip any comments or empty lines right after "Masses"
            while i < len(lines) and (not lines[i].strip() or lines[i].strip().startswith("#")):
                masses_section.append(lines[i])
                i += 1
                
            # Continue reading mass entries until an empty line followed by a non-mass entry
            while i < len(lines):
                line = lines[i]
                if not line.strip():
                    masses_section.append(line)
                    # Check if the next line is the start of a new section
                    if i+1 < len(lines) and lines[i+1].strip() and not lines[i+1].strip()[0].isdigit():
                        break
                elif line.strip() and line.strip()[0].isdigit():
                    # This looks like a mass entry
                    masses_section.append(line)
                else:
                    # This is not a mass entry, we've reached the next section
                    break
                i += 1
            break
        i += 1
    
    # Add all mass lines
    output_lines.extend(masses_section)
    
    # Find and add the Atoms section
    atoms_section = []
    in_atoms = False
    i = 0
    
    while i < len(lines):
        if "Atoms" in lines[i] and "#" in lines[i]:
            in_atoms = True
            atoms_section.append(lines[i])
        elif in_atoms:
            # Stop at the next section or end of file
            if lines[i].strip() and lines[i].strip()[0].isalpha() and lines[i].strip() not in ["Atoms"]:
                if any(section in lines[i] for section in ["Bonds", "Angles", "Dihedrals", "Impropers", "Pair Coeffs", "Bond Coeffs"]):
                    break
            atoms_section.append(lines[i])
        i += 1
    
    # Add all atom lines
    output_lines.extend(atoms_section)
    
    # Write the output file
    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    
    print(f"Processed {input_file} -> {output_file}")

def main():
    # Get current directory
    pristine_dir = os.path.dirname(os.path.abspath(__file__))
    
    if len(sys.argv) < 2:
        print("Usage: python clean_data.py <output_filename> [input_file1] [input_file2] ...")
        print("Example: python clean_data.py clean_structure.data PBI_5_msi2lmp.data")
        print("If no input files are specified, all .data files not starting with 'clean_' will be processed.")
        return
    
    output_filename = sys.argv[1]
    
    if len(sys.argv) > 2:
        # Process specified input files
        for input_filename in sys.argv[2:]:
            input_path = os.path.join(pristine_dir, input_filename)
            if os.path.exists(input_path):
                # Use the provided output filename
                output_path = os.path.join(pristine_dir, output_filename)
                clean_data_file(input_path, output_path)
            else:
                print(f"File not found: {input_path}")
    else:
        # Process all .data files in the directory
        for filename in os.listdir(pristine_dir):
            if filename.endswith('.data') and not filename.startswith('clean_'):
                input_path = os.path.join(pristine_dir, filename)
                # Create output filename by adding 'clean_' prefix
                output_path = os.path.join(pristine_dir, output_filename)
                clean_data_file(input_path, output_path)
                break  # Only process the first file if no specific files provided

if __name__ == "__main__":
    main() 