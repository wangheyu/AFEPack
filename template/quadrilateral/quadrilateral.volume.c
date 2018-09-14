/**我的quadrilateral.tmp_geo中的积分公式的权重加起来不是1， 而是模板单元的面积
 * 而李若的triangle.tmp_geo权重加起来是不是三角形的面积，而是1. 
 * 我能不能通过把这里的面积设置为1.0，而不用修改权重呢。
 */
double quadrilateral_volume(double ** v)
{
	return 1.0; 
};
