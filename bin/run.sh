set -e
export RDWR_DEBUG_VARS="ML_BSSN::gt11 ADMBASE::gxx ML_BSSN::phi"
export RDWR_DEBUG_INDEXES="1 1 1"
./exe/cactus_sim -Roe yes.par
mv CCTK_Proc0.out CCTK_Proc0.yes
./exe/cactus_sim -Roe no.par
mv CCTK_Proc0.out CCTK_Proc0.no
python3 ${0%/*}/trim.py
vimdiff x.yes x.no
