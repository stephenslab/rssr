#ifndef RSSR_TYPES_H
#define RSSR_TYPES_H

#include <RcppEigen.h>

typedef Eigen::Array<double, Eigen::Dynamic, 1> ColumnArray;
typedef Eigen::Array<double, 1,Eigen::Dynamic > RowArray;

typedef Eigen::Ref<const ColumnArray> c_column_internal;
typedef Eigen::Ref<const RowArray> c_row_internal;
typedef Eigen::Ref<const Eigen::ArrayXd> c_arrayxd_internal;
typedef Eigen::Ref<const Eigen::ArrayXi> c_arrayxi_internal;
typedef Eigen::Ref<const Eigen::ArrayXXd> c_arrayxxd_internal;
typedef Eigen::Ref<const Eigen::VectorXd> c_vectorxd_internal;
typedef Eigen::MappedSparseMatrix<double> c_sparseMatrix_internal;

typedef Eigen::Ref<const Eigen::MatrixXd > c_Matrix_internal;
typedef Eigen::Ref<Eigen::ArrayXd> arrayxd_internal;
typedef Eigen::Ref<Eigen::ArrayXXd> arrayxxd_internal;
typedef Eigen::Ref<ColumnArray> column_internal;
typedef Eigen::Ref<RowArray> row_internal;

typedef Eigen::Map<Eigen::ArrayXd> arrayxd_external;
typedef Eigen::Map<Eigen::ArrayXi> arrayxi_external;
typedef Eigen::MappedSparseMatrix<double> sparseMatrix_external;
typedef Eigen::Map<Eigen::MatrixXd > Matrix_external;
typedef Eigen::Map<Eigen::VectorXd> vectorxd_external;

typedef Eigen::Map<Eigen::ArrayXd> mdarray;
typedef Eigen::Map<Eigen::VectorXd> mdvec;
typedef Eigen::Map<const Eigen::VectorXd> c_mdvec;
typedef Eigen::Map<const Eigen::ArrayXd> c_mdarray;

typedef Eigen::Map<Eigen::ArrayXi> miarray;
typedef Eigen::Map<Eigen::ArrayXXd> m2darray;


typedef Eigen::Map<Eigen::MatrixXd> mmat;
typedef Eigen::Map<const Eigen::MatrixXd> c_mmat;



#endif
