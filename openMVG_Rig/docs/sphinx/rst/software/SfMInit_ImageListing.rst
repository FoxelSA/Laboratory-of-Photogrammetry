**********************************
openMVGRig_SfMInit_ImageListing
**********************************

The first task to process an image dataset in openMVG pipelines consist in creating a **sfm_data.json** file that describe the used image dataset.

This structured file lists for each images an object called a View.
This view store image information and lists:

- the image name,
- the image size,
- the internal camera calibration information (intrinsic parameters) (if any)
- the non central camera subposes (R,C) if any.

Each **View** is associated to a camera **Pose** index and an **Intrinsic** camera parameter group (if any). Group means that camera intrinsic parameters can be shared between some Views (that leads to more stable parameters estimation).

.. code-block:: c++

  $ cd build/software/
  $ openMVGRig_SfMInit_ImageListing -i [] -o [] -d [] -m [ ]

Arguments description:

**Required parameters:**

  - **[-i|--imageDirectory]**

  - **[-o|--outputDirectory]**

  - **[-d|--mountPoint]**

  - **[-m|--macAddress]**

**Optional parameters:**





[-b|--upperBound] upper bound limit for timestamp
[-f|--focal] (pixels)
[-g|--gps] GPU/IMU json file

  - **[-c|--channelFile]**  list of subposes to keep for 3D reconstruction

  - **[-r|--rigidRig]** Use non central camera structure (i.e use rigid rig or not)
     -r 0 : no rigid rig
     -r 1 : with rigid rig structure (default)

  - **[-p|--useCalibPrincipalPoint]** Use principal point of calibration or not
     -p 0 : do not use calibration principal point
     -p 1 : use calibration principal point (default)

  - **[-a|--lowerBound]** lower bound limit for timestamp. If image timestamp is less than lower bound, do not use it for 3D reconstruction

  - **[-b|--upperBound]** upper bound limit for timestamp. If image timestamp is bigger than upper bound, do not use it for 3D reconstruction

  - **[-f|--focal]** (value in pixels)

  - **[-g|--gps]** JSON file containing GPS / IMU measurements for each pose.

.. code-block:: c++

  // Example
  $ openMVGRig_SfMInit_ImageListing  -i /home/user/Dataset/ImageDataset_SceauxCastle/images -o /home/user/Dataset/ImageDataset_SceauxCastle/matches -m aa-bb-cc-dd-ee-ff -d /data/

It will produce you a sfm_data.json file that is used by openMVGRig as a scene description.

Once your have computed your dataset description you can compute the image features:

.. toctree::
   :maxdepth: 1

   ./ComputeFeatures.rst

From lists.txt to sfm_data.json
---------------------------------

Old openMVG version (<0.8) use a lists.txt file to describer image parameters.

Example of a lists.txt file where focal is known in advance

.. code-block:: c++

  1404374411_319830-10-RECT-SENSOR.tiff;2592;1936;2059.94;0;1274.91;0;2059.94;967.702;0;0;1;0;10;-0.00667452;0.999943;0.00834353;-0.0144785;-0.00843948;0.99986;0.999873;0.00655278;0.014534;0.0544729;-0.000352938;0.00141571
  1404374411_319830-24-RECT-SENSOR.tiff;2592;1936;2044.67;0;1253;0;2044.67;981.529;0;0;1;0;24;-0.00700488;0.999918;0.0106847;-0.0597896;-0.0110846;0.998149;0.998186;0.00635308;0.0598623;0.0119341;0.80806;0.00110413
  1404374413_319830-10-RECT-SENSOR.tiff;2592;1936;2059.94;0;1274.91;0;2059.94;967.702;0;0;1;1;10;-0.00667452;0.999943;0.00834353;-0.0144785;-0.00843948;0.99986;0.999873;0.00655278;0.014534;0.0544729;-0.000352938;0.00141571
  1404374413_319830-24-RECT-SENSOR.tiff;2592;1936;2044.67;0;1253;0;2044.67;981.529;0;0;1;1;24;-0.00700488;0.999918;0.0106847;-0.0597896;-0.0110846;0.998149;0.998186;0.00635308;0.0598623;0.0119341;0.80806;0.00110413
  ...

You can convert this file to a valid sfm_data.json file by using the **openMVG_main_ConvertList** application.
