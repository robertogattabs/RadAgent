#' Load a set of DICOM studies into a mmButo structure
#' 
#' @description  Allow to load a set of DICOM studies (but just with ONE seriesInstanceUID for each study) into an mmButo structure. 
#' 
#'               N.B: pay attention, using it you c n loo e the sign l 
#' @param obj.mmButo the object mmButo
#' @param pathToOpen the path where the DICOM serie is locate, on the filesystem. It should point to a path where it can find A SET OF FOLDERS, ONE FOR EACH PATIENT
#' @return wel.... nothing. 
#' @examples \dontrun{
#' 
#' obj<-mmButo()
#' MBT.loadCollection(obj.mmButo = obj, pathToOpen='./DICOMSeries/pat001' );
#' 
#' }
#' @export
MBT.loadCollection<-function(obj.mmButo, pathToOpen) {
  obj.mmButo$loadCollection(pathToOpen = pathToOpen )
}

#' Returns the lists of voxels internal to a given ROIName
#' 
#' @description  Once a mmButo object is loaded, with \code{loadCollection} method, it allows to extract all the voxels within a given ROINAME. The behaviour is similar to the \code{GLT.getImageROIVoxels} but is able to work on all the DICOM studies stored in the object.
#' 
#'               N.B: pay attention.... ops! Your wallet is empty! 
#' @param obj.mmButo the object mmButo
#' @param pathToOpen the path where the DICOM serie is locate, on the filesystem. It should point to a path where it can find A SET OF FOLDERS, ONE FOR EACH PATIENT
#' @return wel.... nothing. 
#' @examples \dontrun{
#' 
#' obj<-mmButo()
#' MBT.getROIVoxel(obj.mmButo = obj, ROIName='GTV' );
#' 
#' }
#' @export
MBT.getROIVoxel<-function(obj.mmButo, ROIName) {
  obj.mmButo$getROIVoxel(ROIName = ROIName )
}

#' Retrieve the pixelSpacing values of the stored series
#' 
#' @description  In quickly allows to read the pixel spacing of the series loaded in a mmButo object. Useful when you are in duobt of having differen pixelSpacing and you are not sure you can work without considering the spatial dimension of the voxels..
#' 
#'               N.B: pay attention.... Splat! Ok, bless you!
#' @param obj.mmButo the object mmButo
#' @return a list of arrays; every array contains the x,y,z dimension of the voxel for the interested serie. 
#' @examples \dontrun{
#' 
#' obj<-mmButo()
#' MBT.getROIVoxel(obj.mmButo = obj, ROIName='GTV' );
#' 
#' }
#' @export
MBT.getAllPixelSpacing<-function(obj.mmButo) {
  obj.mmButo$getAllPixelSpacing( )
}

#' Retrieve the imageVoxelCube of the main images
#' 
#' @description  It gives back the complete MRI/CT/PET scans, evantually interpolated to the wished step
#' 
#'               N.B: pay attention.... Splat! Ok, bless you!
#' @param obj.mmButo the object mmButo
#' @param ps.x Optional, default \code{NA}. If you want to interpolate the image to a new step on the x-axes you can specif here the new PixelSpaginc along the x direction
#' @param ps.y Optional, default \code{NA}. If you want to interpolate the image to a new step on the Y-axes you can specif here the new PixelSpaginc along the y direction
#' @param ps.z Optional, default \code{NA}. come on... try to imagine....
#' @return a list of 3D arrays
#' @examples \dontrun{
#' 
#' obj<-mmButo()
#' MBT.getImageVoxel(obj.mmButo = obj);
#' 
#' }
#' @export
MBT.getImageVoxel<-function(obj.mmButo, ps.x = NA, ps.y = NA, ps.z = NA) {
  res<-obj.mmButo$getImageVoxel( ps.x = ps.x, ps.y = ps.y, ps.z = ps.z)
  return(res)
}