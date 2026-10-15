import numpy as np

def read_lammpstrj(filename, atom_type=None, last_n_frames=10, verbose=True):
    frames = []
    frame_count = 0

    with open(filename, "r") as f:
        while True:
            # Read ITEM: TIMESTEP
            line = f.readline()
            if not line:
                break
            if line.strip() != "ITEM: TIMESTEP":
                continue

            timestep = int(f.readline().strip())

            # Read ITEM: NUMBER OF ATOMS
            line = f.readline()  # ITEM: NUMBER OF ATOMS
            natoms = int(f.readline().strip())

            # Read ITEM: BOX BOUNDS
            line = f.readline()  # ITEM: BOX BOUNDS ...
            box = np.zeros((3, 2))
            for dim in range(3):
                parts = f.readline().strip().split()
                box[dim] = [float(parts[0]), float(parts[1])]

            # Read ITEM: ATOMS
            header_line = f.readline().strip()  # ITEM: ATOMS id type x y z ...
            cols = header_line.split()[2:]  # skip "ITEM:" and "ATOMS"

            # Find column indices
            col_map = {name: idx for idx, name in enumerate(cols)}

            has_id = "id" in col_map
            has_type = "type" in col_map
            has_mol = "mol" in col_map
            has_x = "x" in col_map and "y" in col_map and "z" in col_map
            has_xu = "xu" in col_map and "yu" in col_map and "zu" in col_map

            # Read atom data
            atom_data = np.zeros((natoms, len(cols)))
            for i in range(natoms):
                parts = f.readline().strip().split()
                for j, val in enumerate(parts):
                    atom_data[i, j] = float(val)

            # Filter by atom type
            if atom_type is not None and has_type:
                mask = atom_data[:, col_map["type"]].astype(int) == atom_type
                atom_data = atom_data[mask]

            # Extract positions (prefer unwrapped xu,yu,zu if available)
            if has_xu:
                positions = atom_data[:, [col_map["xu"],
                                          col_map["yu"],
                                          col_map["zu"]]]
            elif has_x:
                positions = atom_data[:, [col_map["x"],
                                          col_map["y"],
                                          col_map["z"]]]
            else:
                positions = atom_data[:, :3]

            frame = {
                "timestep": timestep,
                "positions": positions,
                "box": box,
                "natoms": len(positions),
            }

            # Add mol_ids if available
            if has_mol:
                frame["mol_ids"] = atom_data[:, col_map["mol"]].astype(int)

            frames.append(frame)
            frame_count += 1

            if verbose and frame_count % 100 == 0:
                print(f"    Read {frame_count} frames "
                      f"(timestep {timestep})...")

    if verbose:
        print(f"Total frames read: {frame_count}")

    # Keep only last N frames
    if len(frames) > last_n_frames:
        frames = frames[-last_n_frames:]
        if verbose:
            print(f"Keeping last {last_n_frames} frames for analysis.")

    return frames


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        sys.exit(1)

    fname = sys.argv[1]
    atype = int(sys.argv[2]) if len(sys.argv) > 2 else None
    nframes = int(sys.argv[3]) if len(sys.argv) > 3 else 5

    frames = read_lammpstrj(fname, atom_type=atype, last_n_frames=nframes)
    for fr in frames:
        print(f"  Timestep {fr['timestep']}: {fr['natoms']} atoms, "
              f"box z=[{fr['box'][2,0]:.1f}, {fr['box'][2,1]:.1f}]")
