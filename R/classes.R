#	CLASSES.R

setClass("MArrayTC",
   representation("list")
)

setClass("LargeDataObject")
setIs("MArrayTC","LargeDataObject")


