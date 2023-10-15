
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
                
                
def print_diff(old_val: str, new_val: str, linenum: int, 
               col: int, at: int, filename: str) -> None:
    print(f"{filename}: Values do not match at line {linenum}, column {col}:")
    print(f"Current:  {new_val}")
    print((" " * (at + 10))+ ("^" * (len(new_val) - at)))
    print(f"Original: {old_val}\n") 
    

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
                    
        print(f"Detected {diff_count} diffs on file {new_file}")
        input("...Continue [ENTER]\n")
                
            
if __name__ == "__main__":
    SystemExit(compare())
            
            
            
        
                
        
                    
    