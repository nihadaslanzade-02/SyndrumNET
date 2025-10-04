"""
Add new drugs to the pipeline.
"""

import argparse

def main():
    parser = argparse.ArgumentParser(description="Add new drugs")
    parser.add_argument('--lincs', type=str, required=True, help="File with LINCS IDs")
    parser.add_argument('--config', type=str, required=True, help="Config file")
    args = parser.parse_args()
    
    print(f"Adding new drugs from: {args.lincs}")
    print("\nNote: Full implementation requires LINCS L1000 integration.")
    print("This is a placeholder for the extension functionality.")


if __name__ == '__main__':
    main()