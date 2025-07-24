import os
from collections import deque

def extract_last_eigenvalue(filename):
    with open(filename, "rb") as f:
        one_before_last = deque(f, 2)[0]
    
    text = one_before_last.decode("utf-8").strip()
    parts = text.split()
    return parts[-1] if parts else None

print(extract_last_eigenvalue("./eigenvalues.txt"))
