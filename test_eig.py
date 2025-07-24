import os

def extract_last_eigenvalue(filename):
    try:
        with open(filename, "rb") as f:
           try:  # catch OSError in case of a one line file 
              f.seek(-2, os.SEEK_END)
              while f.read(1) != b'\n':
                 f.seek(-2, os.SEEK_CUR)
           except OSError:
              f.seek(0)
           last_line = f.readline().decode() 
           print(last_line)

        # The value is after the last colon or space
        parts = last_line.split()
        return parts
        if parts:
            return float(parts[-1])
    except Exception as e:
        print(f"Failed to read last eigenvalue: {e}")
    return None


print(extract_last_eigenvalue("./eigenvalues.txt"))
