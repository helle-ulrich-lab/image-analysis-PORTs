**Characterization of nuclear foci in *S. cerevisiae* by fluorescence
microscopy**

This repository contains a set of image analysis scripts used in the
study “Processing of DNA polymerase-blocking lesions during genome
replication is spatially and temporally segregated from replication
forks” by Wong et al. (currently under revision at *Molecular Cell*). In
this study, all image analyses were automated using unbiased, customized
scripts written in ImageJ macro and Jython languages with ImageJ FIJI
software (<https://fiji.sc/>). The following analyses were used:

1.  **2D foci counting** (Related to Figures 1D, 1G, S1A and S1B)

    This script determines the number of foci of a fluorescent protein
    in individual nuclei. It analyzes images acquired with three
    channels, but uses only the two fluorescent channels for
    downstream analysis.

    -   Input: Channel 1 (target fluorescent protein), Channel 2 (DAPI) and
        Channel 3 (bright field, not used)

    -   Approach: The DAPI signal (Channel 2) is used to create a
        nuclear mask. Foci (Channel 1) are then segmented in 2D by
        auto-thresholding, and the number of foci per nucleus is determined
        for each nuclear mask.

    -   Critical parameters to optimize: choice of
        auto-thresholding method(s) for nuclei and foci

    -   Output: a table listing the number of foci for each nucleus analyzed

1.  **3D foci counting** (Related to Figures 2A, S2A, S2B and S2C)

    This script determines the number of foci of a fluorescent protein
    in individual nuclei. It performs foci segmentation in 3D to avoid
    fusion of foci in lateral proximity, but in different planes upon
    Z projection. It analyzes images acquired with three channels, but
    uses only the two fluorescent channels for downstream analysis.

    -   Input: Channel 1 (target fluorescent protein), Channel 2 (DAPI) and
        Channel 3 (bright field, not used)

    -   Approach: Nuclear masks are created in 2D from the DAPI signal
        (Channel 2). A “3D TopHat” filter is then applied to the images of
        the foci (Channel 1) before auto-thresholding, and the number of
        foci per nucleus is determined for each nuclear mask.

    -   Critical parameters to optimize: filtering parameter for “3D TopHat”
        filter and choice of auto-thresholding method(s) for nuclei and foci

    -   Output: a table listing the number of foci for each nucleus analyzed

1.  **3D foci counting, colocalization, intensity and volume
    quantification** (Related to Figures 2B, 3B-D, 3G, 4A-D, 5A, 5C-D,
    6C, 6G-H, 6K, S2D-F, S3A-E, S3K, S4A-B, S5A-F, S6B, S6E, S6G)

    This script performs an object-based colocalization analysis for two
    target proteins imaged in two different channels. It determines the
    number of foci per nucleus for both targets, their intensities and
    volumes as well as the degree of overlap with the other target. The
    background fluorescence of one of the targets is used to create
    nuclear masks; if this is not possible, a fluorescent nuclear marker
    should be introduced by other means.

    -   Input: Channel 1 (target fluorescent protein 1 with exclusively
        nuclear localization) and Channel 2 (target fluorescent protein 2)

    -   Approach: Nuclear masks are created from the pan-nuclear background
        signal of target protein 1 (Channel 1), using filtering
        and auto-thresholding. An absolute threshold is applied to each of
        the targets in 3D for foci segmentation, and the number of foci per
        nucleus is determined for each nuclear mask. The co-localization
        value for each object (focus) is calculated by the relative volume
        (in percent) of each object that is covered by the other target. To
        quantify total foci intensity, the 3D foci masks generated from
        segmented foci are used to measure intensities in the original
        fluorescence images (expressed as total integrated densities in
        arbitrary units). Foci volumes are calculated from foci masks
        (expressed as µm^3^ for scaled images or number of voxels for
        unscaled images).

    -   Critical parameters to optimize: minimum size and absolute threshold
        for foci segmentation of each target and choice of auto-thresholding
        method for nuclear segmentation

    -   Output: four tables in total – one table per target protein listing
        the number of foci per nucleus and one table per target protein
        listing the properties of each focus (degree of overlap with the
        other target protein, intensity and volume)

1.  **Zoning assay** (Related to Figure 5B)

    This macro counts foci and determines the location of each focus
    with respect to three nuclear zones of equal volume. It requires a
    fluorescent channel that marks the nuclear volume and a target
    fluorescent protein.

    -   Input: Channel 1 (nuclear signal) and Channel 2 (target
        fluorescent protein)

    -   Approach: 3D nuclear masks are generated by processing the nuclear
        signal (Channel 1) with a “Laplacian of Gaussian” filter
        and auto-thresholding. The resulting 3D objects are tested for
        their circularity. Nuclei that do not pass a threshold value are
        removed before ellipsoid fitting. The fitted ellipsoids are eroded
        into three concentric 3D zones of equal volume. Images of the target
        fluorescent protein (Channel 2) are segmented with an absolute
        threshold, and each object (focus) is assigned to one of the three
        nuclear zones based on the localization of its center.

    -   Critical parameters to optimize: minimum size and absolute threshold
        for foci segmentation of the target fluorescent protein and choice
        of the auto-thresholding method for the nuclear signal

    -   Output: three tables in total – one table listing the number of foci
        per nucleus, one table listing the desired and actual dimensions of
        zones (for quality control), and one table listing the assignment of
        all foci to the nuclear zones

1.  **Quantification of total nuclear intensity** (Related to Figures
    3E, S3F, S3J, S6D)

    This script performs quantification of total and mean nuclear
    fluorescent intensities. It requires a fluorescent channel for
    generating a nuclear mask (DAPI) and a target fluorescent protein.

-   Input: Channel 1 (DAPI) and Channel 2 (target fluorescent protein)

-   Approach: Nuclear masks are created from the DAPI signal (Channel 1)
    using filtering and auto-thresholding. The nuclear masks are then
    used to measure nuclear area (expressed as µm^2^ for scaled image
    and number of pixels for unscaled image) and to quantify total
    intensity of the target fluorescent protein (Channel 2) (expressed
    as total integrated density with an arbitrary unit) and mean
    intensity (total intensity divided by area).

-   Critical parameters to optimize: choice of
    auto-thresholding method(s) for nuclei and foci

-   Output: a table listing total nuclear intensity, area and mean
    intensity per nucleus

1.  **Foci tracking** (Related to Figures 1F, S1D-F)

    This script tracks foci in 2D over time with the plugin TrackMate.

    -   Input: a list of time-lapse movies of a fluorescent protein in 2D
        from individual cells, manually cropped to ensure that there is no
        major shift of cells and the entire nucleus is comprised within the
        z stack

    -   Approach: Foci in a time-lapse movie are identified by a “Laplacian
        of Gaussian” (LOG) detector with filtering and thresholding. Spots
        are connected in time with appropriate constraints. Lifetime
        (in seconds) and number of merging and splitting events are
        determined for each track. Coordinates (expressed as µm for scaled
        image and number of pixels for unscaled image) and mean intensity
        (expressed as an arbitrary unit) are quantified for each focus.

    -   Critical parameters to optimize: parameters for LOG detector
        (expected radius of foci and threshold), constraints for linking and
        splitting events (based on distances and quality of foci between
        frames), and filter for overall track quality

    -   Output: two tables for each movie – one table listing the properties
        of tracks (lifetime, number of merging and splitting events) and one
        table listing coordinates and mean intensities of all foci in each
        frame

