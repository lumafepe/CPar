
FILE_TITLE: str = "cp"

FILENAMES: list[str] = [
    f"{FILE_TITLE}_average.txt",
    f"{FILE_TITLE}_output.txt",
    f"{FILE_TITLE}_traj.xyz"
]

ORIGINAL_FILENAMES: list[str] = [f"output/{name}" for name in FILENAMES]

MATCH: list[tuple[str, str]] = zip(FILENAMES, ORIGINAL_FILENAMES)

 
def is_fp(value_str: str) -> bool:
    
    try: 
        float(value_str)
        return True
    except ValueError as _:
        return False

def at_most_equal(a: str, b: str, level: int = 10) -> bool:
    
    at: int = 0
    
    if len(a) != len(b):
        return False, 0
    
    for char_idx in range(len(a)):
        
        if a[char_idx].isnumeric() and b[char_idx].isnumeric():
            
            if a[char_idx] != b[char_idx] and at <= level:
                return False, at + 1
            
            at += 1
                
    return True, 0
                
def pprint(at: int, string: str, color: str) -> None:
    
    COLOR_CODE: str = "\033[93m" if color == "red" else "\033[92m"
    
    for id, char in enumerate(string):
        if id >= at:
            print(f"{COLOR_CODE}{char}\033[00m", end="")
        else:
            print(char, end="")
    print("")  
               
def print_diff(old_val: str, new_val: str, linenum: int, 
               col: int, at: int, filename: str) -> None:
    
    print(f"{filename}: Values do \033[91m not match \033[00m at line {linenum}, column {col}:")
    print("Current:  ", end="")
    pprint(at, new_val, "red")
        
    print((" " * (at + 10))+ ("^" * (len(new_val) - at)))
    print("Original: ", end="") 
    pprint(at, old_val, "green")
    

def compare():
    
    for (new_file, og_file) in MATCH:
        
        diff_count: int = 0
                
        with open(new_file, "r") as new:
            new_lines: list[str] = new.readlines()[2:]
            
        with open(og_file, "r") as og:
            og_lines: list[str] = og.readlines()[2:]
            
        for i, (new, old) in enumerate(zip(new_lines, og_lines)):
            
            new_values: list[str] = list(filter(lambda v: is_fp(v), new.split(" ")))
            old_values: list[str] = list(filter(lambda v: is_fp(v), old.split(" ")))  
            
            for j, (new_val, old_val) in enumerate(zip(new_values, old_values)):
                cmp, at = at_most_equal(new_val, old_val)
                if not cmp:
                    print_diff(old_val, new_val, i, j, at, new_file)   
                    diff_count += 1
                    
        print(f"Detected \033[91m{diff_count}\033[00m diffs on file {new_file}")
        input("...Continue [ENTER]\n")
                
            
if __name__ == "__main__":
    try: 
        SystemExit(compare())
        
    except KeyboardInterrupt:
        print("")
            
            
            
        
                
        
                    
    
