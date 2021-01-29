# SV utils
Small utilities used in our SV calling pipeline

## FilterSVByValue.py
Annotate SV VCF using any tab-separated file containing var ID and a float value
- Positions of ID column and float value column are specified as arguments
- Float value is annotated in INFO with the specified tag
- A filter tag is added based on float value and max/min specified thresholds

## Fix_svtools_VCF_.py
Fix formatting issues and float GQ values produced by svtool pipeline
- Float GQ values are rounded to int
- Dot values in sample fields are converted to zero
- Poor variants lacking proper FORMAT fields are removed