#! /bin/bash

if [[ $# != 2 ]]; then
  echo "Usage: ${0} dimension prefix"
  echo "  dimension: 维数(1, 2或者3)"
  echo "  prefix: 文件名前缀，包括输入和输出"
  echo "将前面的程序产生的一系列网格文件 prefix*.mesh 转换为 Open DX 的内部格式，并产生 prefix.dx 文件将这些网格文件组合到一起用于显示。"
  exit 1
fi

DIM=$1
PREFIX=$2

echo "object \"group\" class group" >${PREFIX}.dx
value=0
for i in `ls ${PREFIX}?*.mesh`; do 
  dxi=`echo $i | sed -e "s/\.mesh$/.dx/g"`
  mesh2opendx ${DIM} $i ${dxi} 
  echo "  member \"`basename $i`\"" >>${PREFIX}.dx
  echo "  value file \"${dxi}\", \"FEMFunction-${DIM}d\"" >>${PREFIX}.dx
  pnt=`grep "^object 1" ${dxi} | cut -d ' ' -f 12`
  sed -e "s/\(attribute \"ref\" string \"positions\"\)/\1\n\nobject 3 class constantarray type float rank 1 shape 1 item ${pnt} data follows\n${value} \nattribute \"dep\" string \"positions\"/g" ${dxi} | sed -e "s/\(component \"connections\" value 2\)/\1\ncomponent \"data\" value 3/g" >tmp.dx
  mv tmp.dx ${dxi}
  let "value += 1"
done

exit 0

