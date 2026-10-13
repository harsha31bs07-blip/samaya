NAME HAVERLY3
ROWS
 N obj
 E pool_balance
 E pool_quality
 L prod1_quality
 L prod2_quality
 L prod1_demand
 L prod2_demand
COLUMNS
 a obj 6
 a pool_balance 1
 a pool_quality 3
 b obj 13
 b pool_balance 1
 b pool_quality 1
 u obj -9
 u pool_balance -1
 u prod1_quality -2.5
 u prod1_demand 1
 v obj -15
 v pool_balance -1
 v prod2_quality -1.5
 v prod2_demand 1
 z1 obj 1
 z1 prod1_quality -0.5
 z1 prod1_demand 1
 z2 obj -5
 z2 prod2_quality 0.5
 z2 prod2_demand 1
 p obj 0
RHS
 rhs prod1_demand 100
 rhs prod2_demand 200
BOUNDS
 LO bnd p 1
 UP bnd p 3
QCMATRIX pool_quality
 p u -0.5
 u p -0.5
 p v -0.5
 v p -0.5
QCMATRIX prod1_quality
 p u 0.5
 u p 0.5
QCMATRIX prod2_quality
 p v 0.5
 v p 0.5
ENDATA
