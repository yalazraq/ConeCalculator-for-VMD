import mdtraj as md
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import distance

def get_coords_by_resid(traj, frame, res_ids):
    """
    Returns a list of tuples: (ResID, atom_index, [x, y, z]),
    using only the last atom in each residue.
    """
    if frame >= traj.n_frames:
        raise ValueError(f"Frame {frame} is out of range. Max frame index: {traj.n_frames - 1}")

    # Map ResID (resSeq) to atom indices
    resid_to_atom_indices = defaultdict(list)
    for atom in traj.topology.atoms:
        resid = atom.residue.resSeq
        resid_to_atom_indices[resid].append(atom.index)

    result = []
    for resid in res_ids:
        atom_indices = resid_to_atom_indices.get(resid, [])
        if not atom_indices:
            print(f"⚠️ ResID {resid} not found in trajectory!")
            continue

        # Take the last atom in this residue
        last_atom_idx = atom_indices[-1]
        coord = traj.xyz[frame][last_atom_idx] * 10  # Convert from nm to Å
        result.append((resid, last_atom_idx, coord))

    return result

def compute_average_plane(points):
    """Fits an average plane using PCA to approximate the 'floor' of the cone base."""
    centroid = np.mean(points, axis=0)
    cov_matrix = np.cov(points.T)
    _, _, Vt = np.linalg.svd(cov_matrix)
    normal = Vt[2]
    return centroid, normal

def compute_all_triangle_angles(tip, points):
    """Computes all three angles (θ_A, θ_B, θ_C) for each triangle."""
    base_centroid, _ = compute_average_plane(points)
    hypotenuses = [np.linalg.norm(tip - point) for point in points]
    height = min(np.linalg.norm(tip - base_centroid), min(hypotenuses) * 0.9)
    bases = [np.linalg.norm(point - base_centroid) for point in points]
    all_angles = []
    for i, c in enumerate(hypotenuses):
        b = bases[i]
        if c > 0 and height > 0:
            theta_A = np.degrees(np.arcsin(height / c))
            theta_B = np.degrees(np.arcsin(b / c))
            theta_C = 180 - (theta_A + theta_B)
            all_angles.append((theta_A, theta_B, theta_C))
    avg_theta_A = np.mean([angles[0] for angles in all_angles])
    avg_theta_B = np.mean([angles[1] for angles in all_angles])
    avg_theta_C = np.mean([angles[2] for angles in all_angles])
    return all_angles, (avg_theta_A, avg_theta_B, avg_theta_C), base_centroid, height

def plot_final_triangle_with_angles(tip, points, base_centroid, triangle_angles):
    """Visualizes the final triangle structure with all three angles labeled."""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(*tip, color='red', marker='o', s=100, label='Cone Tip (B)')
    ax.scatter(*base_centroid, color='purple', marker='x', s=100, label='Base Centroid')
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], color='blue', label='Base Points (A)')
    for i, point in enumerate(points):
        ax.plot([tip[0], point[0]], [tip[1], point[1]], [tip[2], point[2]], 'k--', label="Hypotenuse (c)" if i == 0 else None)
        ax.plot([tip[0], base_centroid[0]], [tip[1], base_centroid[1]], [tip[2], base_centroid[2]], 'g-', label="Height (h)" if i == 0 else None)
        ax.plot([base_centroid[0], point[0]], [base_centroid[1], point[1]], [base_centroid[2], point[2]], 'b-', label="Base (b)" if i == 0 else None)
        ax.text(point[0], point[1], point[2], f"θ_A={triangle_angles[i][0]:.2f}°", color='blue', fontsize=10)
        ax.text(tip[0], tip[1], tip[2], f"θ_B={triangle_angles[i][1]:.2f}°", color='red', fontsize=10)
        ax.text(base_centroid[0], base_centroid[1], base_centroid[2], f"θ_C={triangle_angles[i][2]:.2f}°", color='green', fontsize=10)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.title("Final Corrected Triangle Visualization with All Angles")
    plt.show()


#Na31-2148H2O.dcd
def runEntireProgram():
    while True:
        try:
            pdbFile = input("Enter Topology file file name (file must be in same directory): ")
            dcdFile = input("Enter DCD file name (file must be in same directory): ")
            traj = md.load_dcd(dcdFile, top=pdbFile)
            print("Total number of frames in current load: " + str(traj.n_frames))
            break
        except Exception as e:
            print("One of the files does not exist or failed to load.")
            print("Error details:", e)
    
    
    all_avg_tip_angles = []
    
    while True:
        base_points = []
        cone_tip = []
        
        while True:
            number = int(input("Enter an integer Frame Number: "))
            frame_index = number-1
            if frame_index < traj.n_frames and frame_index > 0:
                break;
            else:
                print("invalid frame number.\n")
                
    
        # Collect base point coordinates
        while True:
            res_ids_input = []
            try:
                number = int(input("Enter an integer resID for basepoints (or any non-integer to finish): "))
                res_ids_input.append(number)
                coords_info = get_coords_by_resid(traj, frame_index, res_ids_input)
                base_points.append(coords_info[0][2])  
            except ValueError:
                print("Finished collecting basepoints.\n")
                break
            except IndexError:
                print(f"⚠️ ResID {number} not found in trajectory!")
    
        # Collect cone tip coordinate
        tip_ids_input = []
        try:
            while True:
                number2 = int(input("Enter an integer resID for tip: "))
                sureCheck = input("are you sure this is the correct tip resID? (program will crash if resID not in Trajectory). (y/n): ").strip().lower()
                if sureCheck == 'y':
                    break
                    
            tip_ids_input.append(number2)
            tip_info = get_coords_by_resid(traj, frame_index, tip_ids_input)
            cone_tip.append(tip_info[0][2])
        except ValueError:
            print("Invalid input for tip. Skipping this set.\n")
            continue
    
        # Convert to NumPy arrays
        base_points = np.array(base_points)         
        cone_tip = np.array(cone_tip)[0]            
    
        # Compute angles and results
        final_triangle_angles, final_avg_angles, final_base_centroid, final_height = compute_all_triangle_angles(cone_tip, base_points)
    
        # Save average cone tip angle θ_B
        all_avg_tip_angles.append(final_avg_angles[1])
    
        # Print result
        print(f"\n--- Final Average Angle at B ---")
        print(f"Average θ_B (Cone Tip Angle) = {final_avg_angles[1]:.6f}°")
        print(f"--------------------------------\n")
    
        # Plot triangle with angles
        plot_final_triangle_with_angles(cone_tip, base_points, final_base_centroid, final_triangle_angles)
    
        # Ask user if they want to continue
        cont = input("Would you like to analyze another cone? (y/n): ").strip().lower()
        if cont != 'y':
            break
    
    # Final average of all average tip angles
    if all_avg_tip_angles:
        overall_average = sum(all_avg_tip_angles) / len(all_avg_tip_angles)
        print(f"\n=======================================")
        print(f"Final Average of All θ_B Angles: {overall_average:.6f}°")
        print(f"From {len(all_avg_tip_angles)} total iterations.")
        print(f"=======================================\n")
    else:
        print("\nNo angle data was collected.\n")


runEntireProgram()