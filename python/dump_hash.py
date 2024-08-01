#!/usr/bin/env python3

import struct
import sys

def read_uint64_from_file(filename):
    try:
        with open(filename, 'rb') as f:
            #while chunk := f.read(8):  # Read 8 bytes at a time (64 bits)
            chunk = f.read(8)
            while chunk:
                if len(chunk) < 8:
                    print(f"Incomplete data: {chunk}")
                    continue
                value = struct.unpack('<Q', chunk)[0]  # Little-endian unsigned 64-bit integer
                print(value)
                chunk = f.read(8)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 dump_hash.py <binary_file>")
        sys.exit(1)

    binary_file = sys.argv[1]
    read_uint64_from_file(binary_file)

