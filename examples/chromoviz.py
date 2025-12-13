"""
chromoviz.py

Written by Cullen Roth, Ph.D. 

© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Helpers for:
- best-fit rigid alignment between 3D point sets
- 3D spline fitting
- mapping spline evaluation points to genomic coordinates (100 kb bins)
- centromere-based segmentation
- PyVista plotting of two chromosomes side-by-side
"""
## Load in mods 
import numpy as np, pandas as pd, pyvista as pv
from scipy.interpolate import splprep, splev
## Load in glob 
from glob import glob

## Defint sorted glob 
def sortglob(wc):
    return sorted(glob(wc))


def centromere_bins(bin_size: int, centromere,chrom_length: int = None) -> np.ndarray:
    """
    Return the indices of fixed-size bins that intersect a centromere region.

    Parameters
    ----------
    bin_size : int
        Size of each bin in base pairs (e.g., 100_000 for 100 kb).
        Bins are assumed to be:
            bin 0 -> [0, bin_size)
            bin 1 -> [bin_size, 2*bin_size)
            ...
    centromere : (int, int)
        (start, end) genomic coordinates of the centromere in base pairs
        on the same coordinate system as the bins.
        Assumed half-open interval: [start, end)
    chrom_length : int, optional
        Length of the chromosome in base pairs (for sanity checks only).
        If provided, the function clips the centromere interval to [0, chrom_length).

    Returns
    -------
    np.ndarray
        1D array of integer bin indices that intersect the centromere region.
    """

    cstart, cend = centromere

    # Optional: clip centromere to chromosome bounds if chrom_length given
    if chrom_length is not None:
        cstart = max(0, min(cstart, chrom_length))
        cend = max(0, min(cend, chrom_length))
        if cend <= cstart:
            # centromere doesn't overlap chromosome after clipping
            return np.array([], dtype=int)

    # If no chrom_length given, at least ensure a sensible interval
    if cend <= cstart:
        return np.array([], dtype=int)

    # Bin index containing cstart: integer division
    first_bin = cstart // bin_size

    # Bin index containing cend - 1 (last base inside the centromere)
    last_bin = (cend - 1) // bin_size

    # Return all bin indices that intersect [cstart, cend)
    return np.arange(first_bin, last_bin + 1, dtype=int)

# ============================
#  GEOMETRY / ALIGNMENT
# ============================

def best_fit_transform(A: np.ndarray, B: np.ndarray):
    """
    Compute optimal rigid transform (R, t) that maps A -> B in least squares.

    A, B: (N, 3) arrays with point i of A corresponding to point i of B.
    Returns:
        R: (3, 3) rotation matrix
        t: (3,) translation vector
        A_aligned: (N, 3) A after applying R, t
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    assert A.shape == B.shape
    assert A.shape[1] == 3

    centroid_A = A.mean(axis=0)
    centroid_B = B.mean(axis=0)

    AA = A - centroid_A
    BB = B - centroid_B

    H = AA.T @ BB
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Correct for reflection
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = Vt.T @ U.T

    t = centroid_B - R @ centroid_A
    A_aligned = (R @ A.T).T + t

    return R, t, A_aligned


def fit_3d_spline(points, k=3, s=0.0, n_eval=8000):
    """
    Fit a parametric 3D spline through ordered 3D points, and evaluate it.

    points : (N, 3) array
    k      : spline degree (1–5, typically 3)
    s      : smoothing factor (0 = interpolate; >0 = smooth)
    n_eval : number of output spline samples along [0,1]

    Returns:
        spline_points : (n_eval, 3) sampled spline
        u_fine        : (n_eval,) parameter values
        tck           : internal spline representation
    """
    points = np.asarray(points, dtype=float)
    assert points.ndim == 2 and points.shape[1] == 3

    x, y, z = points[:, 0], points[:, 1], points[:, 2]

    # parameterize by cumulative arc length
    d = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2 + np.diff(z) ** 2)
    u = np.concatenate(([0], np.cumsum(d)))
    u /= u[-1] if u[-1] != 0 else 1.0

    tck, _ = splprep([x, y, z], u=u, k=k, s=s)

    u_fine = np.linspace(0, 1, n_eval)
    x_s, y_s, z_s = splev(u_fine, tck)
    spline_points = np.vstack((x_s, y_s, z_s)).T

    return spline_points, u_fine, tck


# ============================
#  GENOMIC COORDINATES
# ============================

def genomic_coords_for_raw_points(num_bins, bin_size_bp=100_000, start_bp=0):
    """
    Genomic midpoint (in bp) for each raw 3D point, assuming each
    point is the midpoint of a fixed-size bin.

    Parameters
    ----------
    num_bins : int
        Number of raw 3D points (bins).
    bin_size_bp : int
        Bin size in base pairs (default 100 kb).
    start_bp : int
        Genomic start position of the first bin (default 0).

    Returns
    -------
    midpoints_bp : (num_bins,) array
        Genomic midpoints in bp for each raw point.
    """
    idx = np.arange(num_bins)
    midpoints_bp = start_bp + (idx + 0.5) * bin_size_bp
    return midpoints_bp


def genomic_coords_for_spline(n_eval, num_bins, bin_size_bp=100_000, start_bp=0):
    """
    Genomic coordinate (in bp) for each spline evaluation point,
    spanning the same genomic interval as the raw 3D bins.

    Parameters
    ----------
    n_eval : int
        Number of spline evaluation points (e.g. 5000).
    num_bins : int
        Number of raw 3D points (bins).
    bin_size_bp : int
        Bin size in base pairs (default 100 kb).
    start_bp : int
        Genomic start position of the first bin (default 0).

    Returns
    -------
    genomic_spline_bp : (n_eval,) array
        Genomic positions in bp for each spline point, from the midpoint
        of the first bin to the midpoint of the last bin.
    """
    chrom_len_bp = num_bins * bin_size_bp
    first_mid_bp = start_bp + 0.5 * bin_size_bp
    last_mid_bp = start_bp + chrom_len_bp - 0.5 * bin_size_bp
    genomic_spline_bp = np.linspace(first_mid_bp, last_mid_bp, n_eval)
    return genomic_spline_bp


# ============================
#  CENTROMERE SEGMENTATION
# ============================

def segment_indices_by_centromere(genomic_pos, cen_start, cen_end):
    """
    Classify each spline index as p-arm, centromere, or q-arm.

    Parameters
    ----------
    genomic_pos : (N,) array
        Genomic coordinate per spline point (bp or Mb or normalized).
    cen_start, cen_end : float
        Centromere interval in same units as genomic_pos.

    Returns
    -------
    p_idx  : array of indices with genomic_pos < cen_start
    cen_idx: array of indices with cen_start <= genomic_pos <= cen_end
    q_idx  : array of indices with genomic_pos > cen_end
    """
    genomic_pos = np.asarray(genomic_pos, float)
    p_mask = genomic_pos < cen_start
    cen_mask = (genomic_pos >= cen_start) & (genomic_pos <= cen_end)
    q_mask = genomic_pos > cen_end

    p_idx = np.where(p_mask)[0]
    cen_idx = np.where(cen_mask)[0]
    q_idx = np.where(q_mask)[0]

    return p_idx, cen_idx, q_idx


# ============================
#  PLOTTING HELPERS
# ============================

def _indices_to_contiguous_segments(idx_array):
    """Convert a sorted array of indices into contiguous [start, end] segments."""
    if len(idx_array) == 0:
        return []
    starts = [idx_array[0]]
    ends = []
    for i in range(1, len(idx_array)):
        if idx_array[i] != idx_array[i - 1] + 1:
            ends.append(idx_array[i - 1])
            starts.append(idx_array[i])
    ends.append(idx_array[-1])
    return list(zip(starts, ends))

# ============================
#  PLOT: COLOR BY CENTROMERE
#  (p-arm / centromere / q-arm)
# ============================

def plot_two_chromosomes_centromere(modelA_spline,modelB_spline,genomic_A,genomic_B,
                                    cenA_start,
                                    cenA_end,
                                    cenB_start,
                                    cenB_end,
                                    tube_radius_A=0.01,
                                    tube_radius_B=0.01,
                                    save_path=None,
                                    window_size=(2400, 1200),
                                    view_vector_A=None,
                                    view_vector_B=None,
                                    zoom_A=1.8,
                                    zoom_B=1.8,
                                    radius_values_A=None,
                                    radius_values_B=None,
                                    n_radius_bins=4,
                                    p_color="darkgrey",
                                    cen_color="white",
                                    q_color="lightblue",
                                    show_centromere=True,
                                    background_color_A="white",
                                    background_color_B="white"
    ):
    """
    Side-by-side plot of two chromosomes colored by p-arm / centromere / q-arm.

    model*_spline : (N,3) spline points
    genomic_*     : (N,) genomic coordinate per spline point (bp or Mb)
    cen*_start/end: centromere interval for that chromosome (same units)
    radius_values_* : optional per-point scores that modulate thickness
    show_centromere : bool If False, do not draw centromere segments (only p- and q-arms).
    """

    #modelA_points = np.asarray(modelA_points)
    modelA_spline = np.asarray(modelA_spline)
    genomic_A = np.asarray(genomic_A, float)

    #modelB_points = np.asarray(modelB_points)
    modelB_spline = np.asarray(modelB_spline)
    genomic_B = np.asarray(genomic_B, float)

    plotter = pv.Plotter(
        shape=(1, 2),
        off_screen=save_path is not None,
        window_size=window_size,
        border=False,
    )

    def add_centromere_colored_spline(
            plotter, spline, genomic_pos,
            cen_start, cen_end,
            base_radius, radius_values=None,
            n_bins=4,
            show_centromere=True, 
            p_color="red", cen_color="yellow", q_color="blue"):

        n = len(spline)
        if genomic_pos.shape[0] != n:
            raise ValueError("genomic_pos length must match spline length")

        p_idx, cen_idx, q_idx = segment_indices_by_centromere(
            genomic_pos, cen_start, cen_end
        )

        if radius_values is not None:
            radius_values = np.asarray(radius_values, float)
            if radius_values.shape[0] != n:
                raise ValueError("radius_values length must match spline length")

            qs = np.linspace(0, 1, n_bins + 1)
            bin_edges = np.quantile(radius_values, qs)
            bin_edges[0] = -np.inf
            bin_edges[-1] = np.inf

            min_r = 0.5 * base_radius
            max_r = 1.5 * base_radius
            bin_radii = np.linspace(min_r, max_r, n_bins)

            def draw_region(idx_array, color):
                if len(idx_array) == 0:
                    return
                vals = radius_values[idx_array]
                bins = np.digitize(vals, bin_edges) - 1
                for b in range(n_bins):
                    mask = bins == b
                    if not np.any(mask):
                        continue
                    group = idx_array[mask]
                    for s, e in _indices_to_contiguous_segments(group):
                        seg_pts = spline[s:e + 1]
                        if len(seg_pts) < 2:
                            continue
                        curve = pv.Spline(seg_pts, len(seg_pts))
                        tube = curve.tube(radius=bin_radii[b])
                        plotter.add_mesh(tube, color=color)

        else:
            def draw_region(idx_array, color):
                if len(idx_array) == 0:
                    return
                for s, e in _indices_to_contiguous_segments(idx_array):
                    seg_pts = spline[s:e + 1]
                    if len(seg_pts) < 2:
                        continue
                    curve = pv.Spline(seg_pts, len(seg_pts))
                    tube = curve.tube(radius=base_radius)
                    plotter.add_mesh(tube, color=color)

        draw_region(p_idx, p_color)
        if show_centromere:            
            draw_region(cen_idx, cen_color)
        draw_region(q_idx, q_color)

    # ----- Chromosome A -----
    plotter.subplot(0, 0)
    plotter.set_background(background_color_A)
    plotter.add_text("Chromosome A", font_size=12)

    #plotter.add_points(
    #    modelA_points,
    #    point_size=point_size,
    #    opacity=opacity,
    #    color="gray",
    #    render_points_as_spheres=True,
    #)

    add_centromere_colored_spline(
        plotter,
        modelA_spline,
        genomic_A,
        cenA_start,
        cenA_end,
        base_radius=tube_radius_A,
        radius_values=radius_values_A,
        n_bins=n_radius_bins,
        p_color=p_color,
        cen_color=cen_color,
        q_color=q_color,
        show_centromere=show_centromere
    )

    plotter.reset_camera()
    if view_vector_A is not None:
        plotter.view_vector(view_vector_A)
    else:
        plotter.view_isometric()
    plotter.camera.zoom(zoom_A)

    # ----- Chromosome B -----
    plotter.subplot(0, 1)
    plotter.set_background(background_color_B)
    plotter.add_text("Chromosome B", font_size=12)

    #plotter.add_points(
    #    modelB_points,
    #    point_size=point_size,
    #    opacity=opacity,
    #    color="gray",
    #    render_points_as_spheres=True,
    #)

    add_centromere_colored_spline(
        plotter,
        modelB_spline,
        genomic_B,
        cenB_start,
        cenB_end,
        base_radius=tube_radius_B,
        radius_values=radius_values_B,
        n_bins=n_radius_bins,
        p_color=p_color,
        cen_color=cen_color,
        q_color=q_color,
        show_centromere=show_centromere
    )

    plotter.reset_camera()
    if view_vector_B is not None:
        plotter.view_vector(view_vector_B)
    else:
        plotter.view_isometric()
    plotter.camera.zoom(zoom_B)

    if save_path is not None:
        plotter.show(screenshot=save_path)
    else:
        plotter.show()

def plot_chromosomes_centromere_grid(
        models_points,
        models_splines,
        genomics,
        cen_starts,
        cen_ends,
        tube_radius=0.01,
        save_path=None,
        window_size=(1024, 1024),
        view_vector=None,
        zoom=1.45,
        radius_values_list=None,
        p_color="tab:green",
        cen_color="w",
        q_color="tab:blue",
        show_centromere=True,
        n_rows=2,
        n_cols=2,
        titles=None,
        background_color="grey",
    ):
    """
    Plot multiple chromosomes in an n_rows x n_cols grid, colored by
    p-arm / centromere / q-arm.

    Parameters
    ----------
    models_points : list of (Ni, 3) arrays
        Noisy point clouds for each chromosome.
    models_splines : list of (Ni, 3) arrays
        Spline points for each chromosome.
    genomics : list of (Ni,) arrays
        Genomic coordinates per spline point (bp or Mb).
    cen_starts, cen_ends : list of float
        Centromere interval for each chromosome.
    tube_radius : float
        Base tube radius.
    radius_values_list : list of (Ni,) arrays or None
        Optional per-point scores for thickness.
    show_centromere : bool
        If False, hide centromere segments.
    n_rows, n_cols : int
        Grid layout.
    titles : list of str or None
        Optional subplot titles.
    """

    models_points = [np.asarray(m) for m in models_points]
    models_splines = [np.asarray(s) for s in models_splines]
    genomics = [np.asarray(g, float) for g in genomics]

    n_chrom = len(models_points)
    if not (len(models_splines) == len(genomics) ==
            len(cen_starts) == len(cen_ends) == n_chrom):
        raise ValueError("All input lists must have the same length.")

    if radius_values_list is None:
        radius_values_list = [None] * n_chrom
    elif len(radius_values_list) != n_chrom:
        raise ValueError("radius_values_list length mismatch.")

    if n_rows * n_cols < n_chrom:
        raise ValueError("Grid too small for number of chromosomes.")

    plotter = pv.Plotter(
        shape=(n_rows, n_cols),
        off_screen=save_path is not None,
        window_size=window_size,
        border=False,
    )

    def draw_chromosome(ax_idx):
        r = ax_idx // n_cols
        c = ax_idx % n_cols
        plotter.subplot(r, c)
        plotter.set_background(background_color)

        pts = models_points[ax_idx]
        spline = models_splines[ax_idx]
        genomic_pos = genomics[ax_idx]
        cen_start = cen_starts[ax_idx]
        cen_end = cen_ends[ax_idx]
        radius_values = radius_values_list[ax_idx]

        #plotter.add_points(
        #    pts,
        #    point_size=point_size,
        #    opacity=opacity,
        #    color="gray",
        #    render_points_as_spheres=True,
        #)

        p_idx, cen_idx, q_idx = segment_indices_by_centromere(
            genomic_pos, cen_start, cen_end
        )

        ## Deal with radi
        if radius_values is not None:
            rv = np.asarray(radius_values, float)
            if rv.shape[0] != spline.shape[0]:
                raise ValueError("radius_values length must match spline length")

            # Turn NaN / ±inf into neutral 0.0
            rv = np.nan_to_num(rv, nan=0.0, posinf=0.0, neginf=0.0)
            finite = np.isfinite(rv)

            if not np.any(finite):
                # No usable values; constant radius everywhere
                radius_per_point = np.full(spline.shape[0], tube_radius, dtype=float)
            else:
                vmax = np.max(np.abs(rv[finite]))
                if vmax == 0:
                    radius_per_point = np.full(spline.shape[0], tube_radius, dtype=float)
                else:
                    # Normalize scores into [-1, 1]
                    vnorm = rv / vmax

                    # Map [-1, 1] -> [min_scale, max_scale] around tube_radius
                    min_scale, max_scale = 0.001, 3
                    scale = (vnorm + 1.0) * 0.5          # -> [0, 1]
                    scale = min_scale + scale * (max_scale - min_scale)
                    radius_per_point = tube_radius * scale
        else:
            # No radius_values: constant radius
            radius_per_point = np.full(spline.shape[0], tube_radius, dtype=float)

        #print(np.min(radius_per_point),np.mean(radius_per_point),np.max(radius_per_point))

        def draw_region(idx_array, color):
            idx_array = np.asarray(idx_array, dtype=int)
            if idx_array.size < 2:
                return

            # Ensure indices are sorted and contiguous subsets are respected
            idx_array = np.sort(idx_array)

            # Extract the subset of spline and radius for this arm
            arm_pts = spline[idx_array]
            arm_radius = radius_per_point[idx_array]

            # Build a polyline for this arm
            curve = pv.Spline(arm_pts, n_points=len(arm_pts))

            # Attach per-point radius as a scalar array
            curve['radius'] = arm_radius.astype(float)

            # Let VTK vary radius by scalar
            tube = curve.tube(radius=tube_radius,scalars='radius')  # radius argument omitted on purpose

            plotter.add_mesh(tube, color=color)

        draw_region(p_idx, p_color)
        if show_centromere:
            draw_region(cen_idx, cen_color)
        draw_region(q_idx, q_color)

        if titles is not None and ax_idx < len(titles):
            if titles[ax_idx] is not None:
                plotter.add_text(titles[ax_idx], font_size=12)

        plotter.reset_camera()
        if view_vector is not None:
            plotter.view_vector(view_vector)
        else:
            plotter.view_isometric()
        plotter.camera.zoom(zoom)

    for i in range(n_chrom):
        draw_chromosome(i)

    if save_path is not None:
        plotter.show(screenshot=save_path)
    else:
        plotter.show() 


def load_fit_align_chromosomes(csv_paths,centromere=(),n_eval=8000,bin_size_bp=100_000,smoothness = 0,mask_list=[],correct_size=0.25):
    """
    Load 3D chromosome models from CSV, fit splines using your existing
    helper functions, align all models to the first, and return:

        - aligned raw points
        - aligned spline points
        - genomic coordinates for spline evaluation points

    Parameters
    ----------
    csv_paths : list of str
        Paths to CSV files. Each must have columns 'x', 'y', 'z'.
    n_eval : int
        Number of evaluation points along each spline (must match what
        your plotting code expects).
    bin_size_bp : int
        Bin size in base pairs for your genomic mapping helper.
    start_bp : int
        Genomic coordinate of the first bin.

    Returns
    -------
    aligned_points_list : list of (Ni, 3) arrays
        Raw points per chromosome, aligned to the first spline.
    aligned_splines_list : list of (n_eval, 3) arrays
        Spline evaluation points, aligned to the first spline.
    genomic_list : list of (n_eval,) arrays
        Genomic coordinates (bp) per spline evaluation point.
        These are the outputs of your genomic mapping helper; they are
        not changed by the rigid-body alignment.
    """

    # ------------------------------------------------------------------
    # 0. Load raw points from the CSVs
    # ------------------------------------------------------------------
    ## Load centromere if mask list was provided
    if len(mask_list):
        assert len(centromere), "[centromere correction] A tuple with interger values must be passed!"
        bins_centormere = centromere_bins(bin_size_bp,centromere)

        correctfrom_bins = [b for b in bins_centormere if b not in mask_list]

    raw_points_list = []
    pos_cols = ["x", "y", "z"]
    for path in csv_paths:
        pts = pd.read_csv(path)[pos_cols]
        ## Edit out points of centroemre or unwanted points
        if len(mask_list):
            ## Mask the bad points in given list
            pts.loc[mask_list,'x'] = pts.loc[max(correctfrom_bins),'x'] + np.random.normal(0,correct_size,len(mask_list))
            pts.loc[mask_list,'y'] = pts.loc[max(correctfrom_bins),'y'] + np.random.normal(0,correct_size,len(mask_list))
            pts.loc[mask_list,'z'] = pts.loc[max(correctfrom_bins),'z'] + np.random.normal(0,correct_size,len(mask_list))

        raw_points_list.append(pts.to_numpy(float))
    ## Set number of bins
    n_bins = len(pts)
    # ------------------------------------------------------------------
    # 1. Fit splines 
    # ------------------------------------------------------------------
    spline_list = []

    for pts in raw_points_list:

        # (1) Fit spline using your spline fitting function
        spline_pts, _, _ = fit_3d_spline(pts,k=3,s=smoothness,n_eval=n_eval)  # <-- REPLACE NAME/SIG
        ## Append to list 
        spline_list.append(spline_pts)

    # ------------------------------------------------------------------
    # 2. Align everything to the first spline using YOUR rigid alignment
    # ------------------------------------------------------------------
    ref_spline = spline_list[0]

    aligned_splines_list = []
    aligned_points_list = []

    for i, (pts, spline_pts) in enumerate(zip(raw_points_list, spline_list)):
        if not i:
            # First is the reference: no transform
            aligned_points_list.append(pts.copy())
            aligned_splines_list.append(ref_spline.copy())
        else:
            ## Align the spline pnts to refernece 
            R, t, spline_aligned = best_fit_transform(spline_pts, ref_spline)
            pts_aligned = (R @ pts.T).T + t

            ## append to lsits 
            aligned_splines_list.append(spline_aligned)
            aligned_points_list.append(pts_aligned)
    
    # ------------------------------------------------------------------
    # 3. Calcualte genomic coordiantes
    # ------------------------------------------------------------------
    ## This is the same repeated list
    genomic_list = [genomic_coords_for_spline(n_eval=ref_spline.shape[0],num_bins=n_bins,bin_size_bp=bin_size_bp) 
                    for i in csv_paths]
    
    ## Return the aligned points 
    return aligned_points_list, aligned_splines_list, genomic_list
