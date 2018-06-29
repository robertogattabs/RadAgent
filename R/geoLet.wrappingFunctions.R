#' Load a DICOM serie into a geoLet object
#' 
#' @description  Allow to load a DICOM serie into a geoLet object.
#' 
#'               N.B: pay attention, it can explode if not properly used
#' @param obj.geoLet the object geoLet to load the DICOM serie in;
#' @param pathToOpen the path where the DICOM serie is locate, on the filesystem
#' @param defaultExtension the 'typical' extension of your DICOM files. Default is '.dcm'
#' @return wel.... nothing. It directly load the serie into the passed geoLet object. Technically it shouldn't work but apparently it does.
#' @examples \dontrun{
#' 
#' obj<-geoLet()
#' GLT.openDICOMFolder(obj = obj, pathToOpen='./DICOMSeries/pat001' );
#' 
#' It is equivalent to:
#' 
#' obj<-geoLet();
#' obj$openDICOMFolder(pathToOpen='./DICOMSeries/pat001' );
#' }
#' @export
GLT.openDICOMFolder<-function(obj.geoLet, pathToOpen, defaultExtension='*.dcm') {
  obj.geoLet$openDICOMFolder(pathToOpen = pathToOpen, defaultExtension=defaultExtension)
}

#' Returns an IMAGE voxel cube
#' 
#' @description  Returns the primary IMAGE voxel cube stored into a geoLet object
#' 
#'               N.B: pay attention, it can cause infertility
#' @param obj.geoLet the object geoLet to load the DICOM serie in;
#' @param ps.x The desided pixelspaxing along the x-axes. If not specifies, by default, it will keep the original pixel spacing. Pay attention: lowest pixel spacings can cause long computation time and crashes due to memory swapping on disk.
#' @param ps.y The desided pixelspaxing along the y-axes. If not specifies, by default, it will keep the original pixel spacing. Pay attention: lowest pixel spacings can cause long computation time and crashes due to memory swapping on disk.
#' @param ps.z come on, genius, try to guess....
#' @return  it returns a 3D array containing ALL the IMAGE voxel cubes. No bullshit.
#' @examples \dontrun{
#' 
#' obj<-geoLet()
#' obj$openDICOMFolder(pathToOpen='./DICOMSeries/pat001' );
#' imageVoxelCube<-GLT.getImageVoxelCube(obj = obj );
#' 
#' }
#' @export
GLT.getImageVoxelCube<-function( obj.geoLet, ps.x=NA, ps.y=NA, ps.z=NA) {
  res<-obj.geoLet$getImageVoxelCube(ps.x = ps.x, ps.y = ps.y, ps.z = ps.z)
  return(res)
}

#' Returns the pixelSpacing along three axes of a geoLet DICOM image serie
#' 
#' @description  Returns the pixelSpacing of a desired DICOM Series previously stored in a geoLet Object
#' 
#'               N.B: pay attention, it can bite
#' @param obj.geoLet the object geoLet to load the DICOM serie in;
#' @param seriesInstanceUID the seriesInstance UID you are interested in. It can also be not specified: in this case it return the pixelSpacing of the 'main' image series (which is normally che interested one)
#' @return  Returns an array containing the pixelspacing along x, y and z
#' @examples \dontrun{
#' 
#' obj<-geoLet()
#' obj$openDICOMFolder(pathToOpen='./DICOMSeries/pat001' );
#' pixelSpacingValue<-GLT.getPixelSpacing(obj = obj );
#' 
#' }
#' @export
GLT.getPixelSpacing<-function( obj.geoLet, seriesInstanceUID = NA) {
  res<-obj.geoLet$getPixelSpacing(seriesInstanceUID = seriesInstanceUID)
  return(res)  
}

#' Returns the list of the ROI available in a geoLet object
#' 
#' @description  Returns the list of the ROI previously loaded in a geoLet object
#' 
#'               N.B: pay attention, it can seduce your wife while you use it
#' @param obj.geoLet the object geoLet to load the DICOM serie in;
#' @return  It returns a matrix where in the second row you can read the ROIName stored into the DICOM RT-struct. Please consider that not alwasy the ROIPointList are really associated to the shown ROIs. This depends from how the DICOM Serie has been exported from TPS.
#' @examples \dontrun{
#' 
#' obj<-geoLet()
#' obj$openDICOMFolder(pathToOpen='./DICOMSeries/pat001' );
#' obj$getROIList( obj.geoLet = obj )
#' 
#' }
#' @export
GLT.getROIList<-function( obj.geoLet ) {
  res<-obj.geoLet$getROIList()
  return(res)  
}

#' Returns the value of a tag of a specified DICOM object
#' 
#' @description  Returns the value of a tag of a specified DICOM object (image file, RTStruct, ...)
#' 
#'               N.B: pay attention, cachexia is one of the most common effect of this function
#' @param obj.geoLet the object geoLet to load the DICOM serie in;
#' @param tag it is a string which indicates the DICOM tag you are interested in. I.e.: "0010,0010" for the PatientID, "0010,0020" for the Patient Name and so on. Please, refer to the DICOM standard for the full list of possible tags. Consider also that this function is designed to return just the most common tags: the most exhostic ones (for example the nested ones) cannot be returned or, anyway, the result is not sure.
#' @param fileName If you know the fileName of the DICOM object you are interested in, you can specify it in this parameter. This parameter is optional: if not specified geoLet will give back, by default, the tag taken from the first image of the main image serie
#'               
#' @return  It returns a matrix where in the second row you can read the ROIName stored into the DICOM RT-struct. Please consider that not alwasy the ROIPointList are really associated to the shown ROIs. This depends from how the DICOM Serie has been exported from TPS.
#' @examples \dontrun{
#' 
#' obj<-geoLet()
#' obj$openDICOMFolder(pathToOpen='./DICOMSeries/pat001' );
#' obj$GLT.getTag( obj.geoLet = obj, tag = "0010,0020")
#' 
#' }
#' @export
GLT.getTag<-function( obj.geoLet , tag ,  fileName="") {
  res<-obj.geoLet$getTag( tag = tag, fileName = fileName)
  return(res)  
}

#' Returns the image voxels internal to a specified ROI
#' 
#' @description  Returns the image voxels internal to a specified ROI. Please consider that this function return the voxel of the main image serie.
#' 
#'               N.B: if you are reading this help probably some thiefs in stoling in your house
#' @param obj.geoLet the object geoLet to load the DICOM serie in
#' @param Structure the ROIName that 'contains' the interested voxels. In order to know which are the structure available, please refers to the \code{GTL.getROIList()} and \code{GTL.getROIPointList} functions.
#' @return  Unfortunately for you, a quite complex structure... It returns a list of tree main elements:
#' \itemize{
#'    \item{ \code{DOM} contains the Dicom Orientation Matrix for each slice of the main image series. director cosines are already appropriately multiplied by pixelspacing value}
#'    \item{ \code{geometricalInformationOfImages} is a list which contains many general information about geometry and space. I.e.: \code{pixelSpacing}, \code{SliceThickness}, \code{SOPClassUID}, \code{Rows}, etc.}
#'    \item{ \code{masked.images} is a list of two elements: \code{voxelCube} and \code{location}. \code{voxelcube} contains the 3D-array of the voxels of the box cropped around the ROI. The value of a cell is NA if the voxel is out of the ROI and the image greylevel value if the voxel is internal to the ROI. Because of the 3D-array is just a 3D box cropped around the ROI, it cannot be located in the original 3D box of the CT scan due to the different dimentions. For this reasons, the second element of the list: \code{location} allow to remap the voxelCube in the biggest voxel cube built on the entire CT scan.  }    
#' }
#' @examples \dontrun{
#' 
#' obj<-geoLet()
#' obj$openDICOMFolder(pathToOpen='./DICOMSeries/pat001' );
#' obj$GLT.getImageROIVoxels( Structure = "GTV" )
#' 
#' }
#' @export
GLT.getImageROIVoxels<-function( obj.geoLet , Structure) {
  res<-obj.geoLet$getROIVoxels( Structure = Structure)
  return(res)  
}

#' Returns the dose voxels internal to a specified ROI
#' 
#' @description  Returns the dose voxels internal to a specified ROI. Please consider that this function return the voxel of the main dose serie. It does not work properly if you have more RD associated to the same SeriesInstanceUID. Because of this function uses mesh calculus, a quite high number of parameters can be passed in order to tune the computation. However, in the most case, the simplest form is preferreable.
#' 
#'               N.B: if you are reading this help probably your car is burning
#' @param obj.geoLet the object geoLet to load the DICOM serie in
#' @param ROIName the ROIName that 'contains' the interested voxels. In order to know which are the structure available, please refers to the \code{GTL.getROIList()} and \code{GTL.getROIPointList} functions.
#' @param newPixelSpacing Optional. By default it uses same pixelValues of the main CT scan: using this parameter you can pass an array (i.e.: \code{c(.9, .9, 1.5 )} ) to interpolate the voxel space by a trilinear interpolation
#' @param plotIT Optional. By default it does not plot anything but if you want, during the comuputation, it can plot the CT scan and the related overlapped dose.
#' @param verbose Optional. \code{FALSE} by default. Set to \code{TRUE} it allow to see some logs during the computation
#' @param forceReCalculus Optional. \code{FALSE} by default. Due to a caching system, implemented to improve performances, If a previous calculus has been interrupted the afterwards computations (on the same ROI)  could have problems. If a computation has been interrupted it is a good practice to set this parameter to \code{TRUE} in the next computation in order to let the algorithm to re-build anything.
#' @param fastEngine \code{TRUE} by default. By default it uses the \code{vcgClostKD()} function which should be quicker but probably a bit less accurate (?). If you want to try something different you can set this parameter on \code{FALSE} and it will use the most classical \code{vcgClost}
#' @param decimation Optional. \code{FALSE} by default. Set it to \code{TRUE} if you want to enable decimation to the mesh structure
#' @param decimationpercentage. Optional. 0.8 by default. If the parameter \code{decimation} is set to \code{TRUE} this parameter indicates the percentage of triangles that should be sacrified in the first (and biggest) mesh-model. Killing triangles allows to have an easier to handle mesh and less memory-consumer but reduce the quality of the approximation. In most cases the optimal value has to be defined empirically. If \code{decimation} is set to \code{FALSE} dont's waste your time in tuning this parameter: decimation will not be performed!
#' @param smoothing Optional, \code{FALSE} by default. This parameter, if se to \code{TRUE} allow to smooth the mesh by a numerical factor indicated in the parameter \code{smoothing.iterations}.
#' @param smoothing.terations Optional, 10by default. Setting this parameter makes sense only if \code{smoothing} was previously set to \code{TRUE}. Increasing the number of iterations we can improve the "smoothing" applied to the original mesh.
#' @param interpolate.dose Optional, \code{TRUE} by default. By default (\code{TRUE}) the returned 3D grid of voxel doses is interpolated according to the pixel spacing of the main CT scan. If you want back a 3D grid voxel doses specified with the same pixelValue but whitout interpolation of the dose value, set this parameter to \code{FALSE}. Personally, I cannot see a good reason to set it to \code{FALSE} but perhaps you need to fit your data with some old algorithm of dose computation and set it to \code{FALSE} can help you in getting more similar results.
#' 
#' @return  It returns a list of tree main elements:
#' \itemize{
#'    \item{ \code{voxelCube.CT} this is a 3D array which contains the voxel cube of the original CT scan }
#'    \item{ \code{voxelCube.Dose} this is a 3D array which contains the voxel cube of the DOSES. The dimension is the same of the \code{voxelCube} because the dose value in the space has been interpolated (via triliear interpolation) in order to overlap the centroids of the voxels between \code{voxelCube.CT} and \code{voxelCube.Dose}. This should allow to quickly plot the expected dose value on a specific pixel/voxel in the CT scan. }
#'    \item{ \code{mesh} is the mesh as returned fro the \pkg{mesh} package.  }    
#' }
#' @examples \dontrun{
#' 
#' obj<-geoLet()
#' obj$openDICOMFolder(pathToOpen='./DICOMSeries/pat001' );
#' obj$GLT.getDoseROIVoxels( Structure = "GTV" )
#' 
#' }
#' @export
GLT.getDoseROIVoxels<-function( obj.geoLet , Structure, 
                                ROIName ,newPixelSpacing=NA, plotIT = FALSE, 
                                verbose=FALSE, forceReCalculus=FALSE, fastEngine = TRUE,
                                decimation=FALSE, decimation.percentage=0.8, 
                                smoothing=FALSE, smoothing.iterations = 10,
                                interpolate.dose = TRUE ) {
  res<-obj.geoLet$extractDoseVoxels( Structure = Structure,
                                     ROIName ,newPixelSpacing=newPixelSpacing, plotIT = plotIT, 
                                     verbose=verbose, forceReCalculus=forceReCalculus, fastEngine = fastEngine,
                                     decimation=decimation, decimation.percentage=decimation.percentage, 
                                     smoothing=smoothing, smoothing.iterations = smoothing.iterations,
                                     interpolate.dose = interpolate.dose   )
  return(res)  
}

#' Returns the DVH calculated from data stored into a geoLet object
#' 
#' @description  This function calculates DVH starting from DICOM RT-DOSE objects, previously loaded into a geoLet object.
#' 
#'               N.B: if you are reading this help, your life is going to end in 10, 9, 8, 7, ...
#' @param obj.geoLet the object geoLet to load the DICOM serie in
#' @param ROIName the name of the ROI you are interested in having the DVH
#' @param newPixelSpacing Optional. By default it uses same pixelValues of the main CT scan but if you want to interpolate the space in order to have a bigger detail you can increase the precision. Pay attention: it can take a lot o memory!!!
#' @param justTheDVH Optional. \code{TRUE} by default. Because of this function is also able to calculate other ancillary data structures, if you want you can set this parameter to \code{FALSE} and get also the IMAGE voxel cube and te DOSE voxel cube (sampled at the same points in the 3D space)
#' @param verbose Optional. \code{FALSE} by default. Set to \code{TRUE} it allow to see some logs during the computation
#' @param forceReCalculus Optional. \code{FALSE} by default. Due to a caching system, implemented to improve performances, If a previous calculus has been interrupted the afterwards computations (on the same ROI)  could have problems. If a computation has been interrupted it is a good practice to set this parameter to \code{TRUE} in the next computation in order to let the algorithm to re-build anything.
#' @param fastEngine \code{TRUE} by default. By default it uses the \code{vcgClostKD()} function which should be quicker but probably a bit less accurate (?). If you want to try something different you can set this parameter on \code{FALSE} and it will use the most classical \code{vcgClost}
#' @param decimation Optional. \code{FALSE} by default. Set it to \code{TRUE} if you want to enable decimation to the mesh structure
#' @param decimationpercentage. Optional. 0.8 by default. If the parameter \code{decimation} is set to \code{TRUE} this parameter indicates the percentage of triangles that should be sacrified in the first (and biggest) mesh-model. Killing triangles allows to have an easier to handle mesh and less memory-consumer but reduce the quality of the approximation. In most cases the optimal value has to be defined empirically. If \code{decimation} is set to \code{FALSE} dont's waste your time in tuning this parameter: decimation will not be performed!
#' @param smoothing Optional, \code{FALSE} by default. This parameter, if se to \code{TRUE} allow to smooth the mesh by a numerical factor indicated in the parameter \code{smoothing.iterations}.
#' @param smoothing.terations Optional, 10by default. Setting this parameter makes sense only if \code{smoothing} was previously set to \code{TRUE}. Increasing the number of iterations we can improve the "smoothing" applied to the original mesh.
#' 
#' @return  Depending on the status of the parameter \code{justTheDVH}. If it is set to \code{TRUE} the function returns a \code{dvhmatrix} object. If set to \code{FALSE} it returns a list of three elements:
#' \itemize{
#'    \item{ \code{voxelCube.CT} this is a 3D array which contains the voxel cube of the original CT scan }
#'    \item{ \code{voxelCube.Dose} this is a 3D array which contains the voxel cube of the DOSES. The dimension is the same of the \code{voxelCube} because the dose value in the space has been interpolated (via triliear interpolation) in order to overlap the centroids of the voxels between \code{voxelCube.CT} and \code{voxelCube.Dose}. This should allow to quickly plot the expected dose value on a specific pixel/voxel in the CT scan. }
#'    \item{ \code{DVHobj}  a \code{dvhmatrix} object }    
#' }
#' @examples \dontrun{
#' 
#' obj<-geoLet()
#' obj$openDICOMFolder(pathToOpen='./DICOMSeries/pat001' );
#' a<-obj$GLT.calculateDVH( Structure = "GTV" )
#' 
#' }
#' @export
GLT.calculateDVHs<-function( obj.geoLet , ROIName,
                             newPixelSpacing=NA , justTheDVH=TRUE,
                             verbose=FALSE, forceReCalculus=FALSE, fastEngine = TRUE,
                             decimation=FALSE, decimation.percentage=0.8, 
                             smoothing=FALSE, smoothing.iterations = 10  ) {
  res<-obj.geoLet$calculateDVH( ROIName = ROIName,
                                newPixelSpacing=newPixelSpacing , justTheDVH=justTheDVH,
                                verbose=verbose, forceReCalculus=forceReCalculus, fastEngine = fastEngine,
                                decimation=decimation, decimation.percentage=decimation.percentage, 
                                smoothing=smoothing, smoothing.iterations = smoothing.iterations  )
  return(res)  
}

#' Remap the ROI points according to the raster geomery of image
#' 
#' @description  Remap the 3D points of the chosen ROI in the 3D bitmap geometry of the main image
#' 
#'               N.B: if you are reading this help, Clostridium difficile will bless you...
#' @param obj.geoLet the object geoLet to load the DICOM serie in
#' @param ROIName the name of the ROI you are interested in remapping
#' @return It return a list of points located in the bitmap geometry of the CT (or MRI) scan. By this function you can in a quite easy way identify ROI points on the image slices
#' @examples \dontrun{
#' 
#' obj<-geoLet()
#' obj$openDICOMFolder(pathToOpen='./DICOMSeries/pat001' );
#' a<-obj$GLT.calculateDVH( ROIName = "GTV" )
#' 
#' }
#' @export
GLT.rotateToAlign<-function( obj.geoLet , ROIName) {
  res<-obj.geoLet$rotateToAlign( ROIName = ROIName )
  return(res)  
}