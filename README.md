# MainBeltComets

This repository was built to detect active objects in the main asteroid belt from images in the OSSOS survey.
The majority of scripts with this purpose are contained in the folder 'getImages'

In 'find_family.py' a list of asteroids can be generated based on status (i.e. 0-3 for not in family, inner, main and outer belt) or by a specific family designation( e.g. all object in the Hungaria family '434'). The list is created based on parsing the AstDys database. The output can be found in the folder getImages/family_lists.

In 'get_images.py' a SSOIS query is made for each object in the asteroid list as output from the script above. This query is for the images from the OSSOS survey. A list of objects with exposures in this data set is output along with information about the coordinates of the object, the exposure time, etc. The output can be found in the folder getImages/image_lists.

In 'get_stamps.py' the list of exposures is used to cut postage stamps of the asteroid from each exposure and upload them to VOSpace in the folder vos:kawebb/postage_stamps/all if the asteroid has a family designation or vos:kawebb/postage_stamps/none if it does not. This organization was chosen arbitrarily.

In 'sep_phot.py' the postage_stamps are run through photometry software (SEP) and the precise asteroid location is determined from its position relative to the predicted location, degree of trailing, and magnitude. Information regarding the photometry and coordinates if output as familyname_all_output.txt and similar but less information is output as familyname_output.txt. The reason two files are output is that not all the information in the first file is relevant for continuing along in the study, but may be useful for ensuring the pipeline is working properly. The output can be found in the folder getImages/phot_output.
To run this for a list of asteroids with output in the format as written by 'get_images.py', the script 'do_all.py' can be used. Running 'sep_phot.py' outright is typically best for looking at a particular asteroid and exposure. 

In 'mbc_detection.py' a PSF is generated from the asteroid stamp and compared to a stellar model PSF generated from the OSSOS MOP. As of yet this comparison is not functional.

In  'parse_for_object_metadata.py' the orbital information of each asteroid that we have an exposure of is collected from the AstDys database. This is used for plotting useful to see which orbital phase space we can detect the most asteroids in.

In 'parse_is_output.py' the photometry output is used to generate light curves. As of yet this is not fully functional.