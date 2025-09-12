Fibro=["OU-HiC-C2-pfib-D16","OU-HiC-C2-pfib-D4","OU-HiC-C3-pfib-D16","OU-HiC-C3-pfib-D4",
       "OU-HiC-N6-pfib-D16","OU-HiC-N6-pfib-D4","OU-HiC-N9-pfib-D16","OU-HiC-N9-pfib-D4"]

IPSC=["OU-HiC-C2IPSC1-R1","OU-HiC-C2IPSC2-R1","OU-HiC-C3IPSC3-R2","OU-HiC-C3IPSC4-R1",
      "OU-HiC-N6IPSC1-R1","OU-HiC-N6IPSC2-R1","OU-HiC-N9IPSC1-R2","OU-HiC-N9IPSC2-R1"]

TimeCourse=["OU-HiC-N6IPSC2-EP-W0","OU-HiC-N6IPSC2-MN-W0","OU-HiC-N6IPSC2-MN-W3","OU-HiC-N6IPSC2-MN-W6",
            "OU-HiC-N9IPSC2-EP-W0","OU-HiC-N9IPSC2-MN-W0","OU-HiC-N9IPSC2-MN-W3","OU-HiC-N9IPSC2-MN-W6"]

MN=["OU-HiC-C2IPSC1MND45-R1","OU-HiC-C2IPSC2MND45-R1","OU-HiC-C3IPSC3MND45-R1","OU-HiC-C3IPSC4MND45-R1",
    "OU-HiC-N6IPSC1MND45-R2","OU-HiC-N6IPSC2MND45-R2","OU-HiC-N9IPSC1MND45-R2","OU-HiC-N9IPSC2MND45-R1",
    "OU-HiC-N9IPSC2MND45-R2"]
indi=Fibro+IPSC+TimeCourse+MN

Fibro_G=["OU-HiC-ALS-Fib-All","OU-HiC-ALS-Fib-D16","OU-HiC-ALS-Fib-D4",
         "OU-HiC-Control-Fib-All","OU-HiC-Control-Fib-D16","OU-HiC-Control-Fib-D4"]
IPSC_G=["OU-HiC-ALS-IPSC-All","OU-HiC-C2-IPSC","OU-HiC-C3-IPSC",
        "OU-HiC-Control-IPSC-All","OU-HiC-N6-IPSC","OU-HiC-N9-IPSC"]
TimeCourse_G=["OU-HiC-N9-N6-EP-W0","OU-HiC-N9-N6-MN-W0",
              "OU-HiC-N9-N6-MN-W3","OU-HiC-N9-N6-MN-W6"]
MN_G=["OU-HiC-ALS-iPSN-D45","OU-HiC-C9-2-iPSN-D45","OU-HiC-C9-3-iPSN-D45",
    "OU-HiC-Control-iPSN-D45","OU-HiC-N6-iPSN-D45","OU-HiC-N9-iPSN-D45"]
grps=Fibro_G+IPSC_G+TimeCourse_G+MN_G

timeColors = {
    'OU-HiC-N6IPSC2_EP_W0' : '#800002',
    'OU-HiC-N6IPSC2_MN_W0' : '#6ED2FF',
    'OU-HiC-N6IPSC2_MN_W3' : '#0080FF',
    'OU-HiC-N6IPSC2_MN_W6' : '#0000FF',
    'OU-HiC-N9IPSC2_EP_W0' : '#800002',
    'OU-HiC-N9IPSC2_MN_W0' : '#6ED2FF',
    'OU-HiC-N9IPSC2_MN_W3' : '#0080FF',
    'OU-HiC-N9IPSC2_MN_W6' : '#0000FF',
}

chrms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
'chr21','chr22','chrX','chrY','chrM']














