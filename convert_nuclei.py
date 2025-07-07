import re

def parse_cpp_shell_states(text_definition: str):
    """Parses shell states from a C++ initializer list like '{{n,l},{n,l}}'."""
    # Regex to find all occurrences of {d,d}
    matches = re.findall(r'\{(\d+),(\d+)\}', text_definition)
    return [(int(match[0]), int(match[1])) for match in matches]

def format_to_shellstate_list(particle_type: str, states: list, items_per_line: int = 8):
    """Formats a list of states into the target ShellState format."""
    
    if not states:
        return f"                {particle_type} = [ & \n                ]"

    # Create all the "ShellState(n,l)" strings
    shell_strings = [f"ShellState({s[0]},{s[1]})" for s in states]

    output_lines = [f"                {particle_type} = [ &"]
    
    # Group the strings into lines
    for i in range(0, len(shell_strings), items_per_line):
        chunk = shell_strings[i:i+items_per_line]
        
        # Split the chunk into two halves for formatting
        midpoint = len(chunk) // 2 if len(chunk) % 2 == 0 else (len(chunk) + 1) // 2
        part1 = ", ".join(chunk[:midpoint])
        part2 = ", ".join(chunk[midpoint:])
        
        line = "                    " + part1
        if part2:
            line += ", " + part2
        
        output_lines.append(line + ", &")

    # Clean up the last line to match the requested format
    if len(output_lines) > 1:
        # Remove the final comma and trailing space from the last content line
        last_line = output_lines[-1]
        output_lines[-1] = last_line.rsplit(',', 1)[0] + " &"

    output_lines.append("                ]")
    return "\n".join(output_lines)


def translate_nuclei_in_cpp(cpp_content: str):
    """
    Finds all nucleus definitions in the C++ source and translates them
    to the ShellState format.
    """
    
    full_output = ""
    
    # Regex to find each if(nuclType==...) block
    nucleus_blocks = re.findall(r'if\(nuclType==(\d+)\)\{([\s\S]*?)\n    \}', cpp_content)
    
    for nucl_type, block_content in nucleus_blocks:
        full_output += f"// Converted format for Nucleus Type: {nucl_type}\n"
        
        # Find the text definition for protons
        protons_match = re.search(r'protons\s*=\s*(\{[\s\S]*?\});', block_content)
        
        proton_states = []
        if protons_match:
            proton_text_def = protons_match.group(1)
            proton_states = parse_cpp_shell_states(proton_text_def)
            full_output += format_to_shellstate_list("protons", proton_states) + "\n"

        # Find the text definition for neutrons - improved regex to handle multi-line definitions
        neutrons_match = re.search(r'neutrons\s*=\s*([\s\S]*?);', block_content)
        
        if neutrons_match:
            neutron_def = neutrons_match.group(1).strip()
            # Handle the case where neutrons are defined as a copy of protons
            if neutron_def == 'protons':
                neutron_states = proton_states
            else:
                neutron_states = parse_cpp_shell_states(neutron_def)
            
            full_output += format_to_shellstate_list("neutrons", neutron_states) + "\n\n"
            
    return full_output.strip()

if __name__ == "__main__":
    try:
        with open('true.cpp', 'r') as f:
            cpp_file_content = f.read()
        
        translated_output = translate_nuclei_in_cpp(cpp_file_content)
        
        print(translated_output)
        
    except FileNotFoundError:
        print("Error: true.cpp not found. Make sure it's in the same directory as this script.")