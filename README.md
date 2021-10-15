# 3D+t Primordia Segmentation and Tracking

## Preprocessing
- `ConvertH5ToTiff.m`: Expects an `*.h5` file as input and extracts `*.tif` files for each channel and time point. This is required for subsequent processing of the individual images. Converting an `*.nd2` file to `*.h5` can be performed with Fiji (www.fiji.sc).

- `PerformImageNormalization.m`: Use this script to perform a quantile-based intensity normalization of the images. Expects an input folder containing individual `*.tif` files for each frame and an output folder to write the result images to. If needed, the quantile to use for normalization can be adjusted in line 37 of the script.

- `PerformElastixRegistration.m`: This script allows you to precompute a deformable registration from frame to frame that can assist the tracking algorithm. The step is optional and the remaining scripts also work without the precomputed registration. If used, provide the folder containing the raw nucleus images as an input and a directory to save the output transformations to. If elastix is not found, make sure to install it and to double-check the elastix path in line 34 of the script. Elastic is open-source and can be obtained from https://elastix.lumc.nl/.

## Data Synthesis (Optional)

This step is only required if no pretrained models are available yet. The idea is to create synthetic training images based on a small database of template nuclei that are altered in a randomized fashion to mimick different shapes and intensity appearance that would be expected to be observed in real data.

1. First, create a few 3D images with segmentations of the nuclei. If you want to do it the hard way, manually create some segmentations or alternatively, just use an automatic method like TWANG that potentially contains a few inaccuracies but mostly is sufficient for generating the training data segmentation.

2. Use the segmented images to create a training image data base with the script `CreateSyntheticNucleiImageDatabase.m`. This script requests the input folder for the segmentation images and additionally the folder with the corresponding raw images of the nuclei. Make sure that both folders contain the same number of images with a naming scheme that preserves the ordering of the files. This is crucial to have matching pairs of segmentation and raw images. Parameters that can optionally be adjusted are the background intensity threshold (`intensityThreshold`), the window size used for assessing the region statistics (`windowSize`) and the padding used for cropping nucleus candidates from the images (`regionPadding`).

3. Once the database is created it can be used to generate new images with similar shape and intensity properties. Use the script `CreateSyntheticNucleiImageFromDatabase.m` to accomplish this. There is a bunch of parameters at the beginning of the script that can be used to modify the data generation process.

4. Use the generated images to train a Cellpose3D model (https://github.com/stegmaierj/Cellpose3D).

## Segmentation

1. Use a pretrained Cellpose3D model to generate the flow fields for each of the nucleus input images. The flow fields should be stored as separate files, i.e., for each raw nucleus input image it is expected to have an additional 3D image fro the x, y, z gradients respectively.

2. Use XPIWIT (https://bitbucket.org/jstegmaier/xpiwit/downloads/) for performing the segmentation on the basis of a simple binary threshold that separates foreground and background and using the flow fields to perform the splitting of the nuclei. Make sure to use one of the latest versions, as the Cellpose3D support was only added recently. Load the XPIWIT pipeline `OtsuThresholdAndCellpose_PrimordiumCellSegmentation.sav`, specify the inputs in the XPIWITGUI and process the images to obtain the nucleus channel instance segmentation.

## Segmenation Procedure

1. Compute deformations of the images in a backwards fashion, i.e., setting t-1 as the fixed image and t as the moving image. Can be performed automatically using elastix as implemented in the script PerformElastixRegistration.m

2. Perform initial detection and segmentation on the original nuclei and membrane images. Use the XPIWIT pipelines NucleiEnhancementDeploy.xml (includes Otsu based foreground detection) or NucleiEnhancedDeployWithoutOtsu.xml (does not filter the detections, i.e., more reliable w.r.t. intensity variations)

3. Perform the initial tracking based on the identified transformations from the registration as well as the detection/segmentation results. Initially, everything is tracked and later the masks of the maximum intensity projections are used to filter the desired tracks. Can be done with the script PerformIterativeSeedBasedNNTracking.m or PerformIterativeSegmentationBasedTracking.m to either perform the tracking using hierarchical clustering or using an extended NN method. Both methods operate in a backwards fashion.

4. Open the initial raw track project trackedProject.prjz in SciXMiner and select only those tracks that intersect with the binarized maximum intensity projection in XY and XZ direction using the script "2_ReidentifyMeanIntensity.m".

5. The newly generated output variable is then used by the script "3_CreateCleanedTrackedImages.m" to generate the nuclei segmentation images only for the selected tracks.

6. Finally, the tracked nuclei images are used to propagate the labels of the nuclei to the membrane segmentation that is obtained with Dennis' Cellpose approach. This can be done with the script "4_ConvertSegmentationTiffsToTracked.m".

## Tracking


## Track Correction
1. Create BigDataViewer representation of the image data. Open up the raw image sequence in Fiji as a hyperstack. If using individual `*.tif` files, just drag and drop the folder containing all frame images to Fiji and import the entire series. To convert the currently open hyperstack to the BigDataViewer (BDV) format, use the BigDataViewer plugin from `Plugins -> BigDataViewer -> Export Current Image as XML/HDF5` and remember the save location. The XML file  will later be needed to setup the Mastodon project.

2. Convert the tracked tiff images to TGMM format that can in turn be imported to Mastodon. To do this, simply run the script `ConvertTrackedImagesToTGMM.m`. This will first ask for an output folder where the TGMM tracking files should be saved to and the for an input folder. The input folder is expected to contain individual 3D stacks for each time point, where the contents of the 3D stacks are label images that are tracked over time (e.g. the images generated in the tracking step described above). TGMM is only used as an auxiliary format here and will be converted to Mastodon, which actually serves as tool for the manual correction of the tracking data.

3. Create a new Mastodon project by first starting Mastodon using `Plugins -> Mastodon`. Next, create new Mastodon project using the TGMM importer functionality `import TGMM`. For the entry `Browse to a BDV file (xml/h5 pair)` select the raw data that you converted to the BDV format in step 1. Specify the TGMM folder that contains your exported XML files. Adjust the timepoint pattern to match your file names. For instance, if the files are named according to the following scheme `20210506_Stiching_EGFP_Caax-H2a-mCherry_t=003_RegionProps.xml`, search for the part of the filename that contains the time code and replace the frame number `003` by a `%03d`. The importer will then try to load the files according to the provided scheme or complain if the files are not found.

4. If the import was successful, you can start editing the tracks as needed using the standard Mastodon functionality. See here for a detailed documentation on how to use the tool https://github.com/mastodon-sc/mastodon . It is sufficient to just focus on the cells that you actually want to include in your final analysis, i.e., not all of the tracks need to be corrected. Primarily focus on undersegmentation errors by deleting the existing undersegmented node, by adding two nodes at the correct locations and then linking the fragmented tracks properly.

5. Once done with the track editing, simply save the files to (MaMuT/Mastodon) and use the script `ConvertMaMuTAnnotationsToImage.m`. The script will ask for an image input folder, a MaMuT project file (`*_mamut.xml`) and an output folder to save the result images to. The tracking IDs of the Mastodon track correction will then be propagated to the image data and the features for all contained cells are directly exported as `*.csv` files as well. This combination of 3D label images and `*.csv` files will then serve as input for the final analysis steps.