#!/bin/bash

DATADIR="../../data/small-all"
LOGDIR="small-all"
rm -f parallelcommands

# SOLVER="scip"
SOLVER="gurobi"
JOURNAL=neurips2022

for j in $DATADIR/*.dat
do
    echo $j

    # # ALL 3
    # export MODELNAME="ssn2n1nI"
    # export MODEL="smallSets.mod v.mod compact.mod n2.mod n1nI.mod end.mod"
    # OUTPUT=../../log/$JOURNAL/$LOGDIR/`basename $j`-$MODELNAME-$SOLVER.out
    # echo "ampl $MODEL $j $SOLVER.run > $OUTPUT" >> parallelcommands

    # # 2 + 1
    # export MODELNAME="ssn2n1"
    # export MODEL="smallSets.mod vAntisymm.mod compact.mod n1.mod n2.mod end.mod"
    # OUTPUT=../../log/$JOURNAL/$LOGDIR/`basename $j`-$MODELNAME-$SOLVER.out
    # echo "ampl $MODEL $j $SOLVER.run > $OUTPUT" >> parallelcommands

    # # 2 + inf
    export MODELNAME="ssn2nI"
    export MODEL="smallSets.mod v.mod compact.mod n2.mod nI.mod end.mod"
    OUTPUT=../../log/$JOURNAL/$LOGDIR/`basename $j`-$MODELNAME-$SOLVER.out
    echo "ampl $MODEL $j $SOLVER.run > $OUTPUT" >> parallelcommands

    # 2
    # export MODELNAME="ssn2"
    # export MODEL="smallSets.mod vAntisymm.mod compact.mod n2.mod end.mod"
    # OUTPUT=../../log/$JOURNAL/$LOGDIR/`basename $j`-$MODELNAME-$SOLVER.out
    # echo "ampl $MODEL $j $SOLVER.run > $OUTPUT" >> parallelcommands

    # # 1+inf
    # export MODELNAME="ssn1nI"
    # export MODEL="smallSets.mod v.mod compact.mod n1nI.mod end.mod"
    # OUTPUT=../../log/$JOURNAL/$LOGDIR/`basename $j`-$MODELNAME-$SOLVER.out
    # echo "ampl $MODEL $j $SOLVER.run > $OUTPUT" >> parallelcommands

    # # inf
    # export MODELNAME="ssnI"
    # export MODEL="smallSets.mod v.mod compact.mod nI.mod end.mod"
    # OUTPUT=../../log/$JOURNAL/$LOGDIR/`basename $j`-$MODELNAME-$SOLVER.out
    # echo "ampl $MODEL $j $SOLVER.run > $OUTPUT" >> parallelcommands

    # # 1
    # export MODELNAME="ssn1"
    # export MODEL="smallSets.mod vAntisymm.mod compact.mod n1.mod end.mod"
    # OUTPUT=../../log/$JOURNAL/$LOGDIR/`basename $j`-$MODELNAME-$SOLVER.out
    # echo "ampl $MODEL $j $SOLVER.run > $OUTPUT" >> parallelcommands

done


# ./run.sh $DATADIR $SOLVER >> parallelcommands

# export MODELNAME="ssn2n1"
# export MODEL="smallSets.mod vAntisymm.mod compact.mod n2.mod n1.mod end.mod"
# ./run.sh $DATADIR $SOLVER >> parallelcommands

# export MODELNAME="ssn2n1nI"
# export MODEL="smallSets.mod v.mod compact.mod n2.mod n1nI.mod end.mod"
# ./run.sh $DATADIR $SOLVER >> parallelcommands

# run as the last one... expecting it to halt at timeout for almost all the instances
# ./run.sh $DATADIR $SOLVER >> parallelcommands

# export MODELNAME="fsn2nI"
# export MODEL="fullSets.mod v.mod compact.mod n2.mod nI.mod end.mod"
# ./run.sh $DATADIR $SOLVER >> parallelcommands

# export MODELNAME="ssn2nIpo"
# export MODEL="smallSets.mod v.mod compact.mod n2.mod nI.mod packingOrbitope.mod end.mod"
# ./run.sh $DATADIR $SOLVER >> parallelcommands

# SOLVER="cplex"



