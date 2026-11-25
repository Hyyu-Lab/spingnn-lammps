"""Independent-chain GPST inference from stored sufficient moments."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent / 'reference'))
import tspin


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records', type=Path, nargs='+', required=True)
    parser.add_argument('--seed-ids', nargs='+', required=True)
    parser.add_argument('--kb', type=float, required=True, help='Native kB in eV/K')
    parser.add_argument('--dt', type=float, required=True, help='Saved-frame interval in ps')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    records = []
    for path in args.records:
        with np.load(path, allow_pickle=False) as data:
            records.append({key: data[key] for key in data.files})
    result = tspin.analyze(records, args.seed_ids, kb=args.kb, dt=args.dt)
    def serial(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)
    args.output.write_text(json.dumps(result, default=serial, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
