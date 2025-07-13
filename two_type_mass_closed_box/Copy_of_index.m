
gb=1e-3;

fparam=1e2;

ic1.m=1;
ic1.gb=gb;

ic2.m=1836;
ic2.gb=gb;

N=1e6;
max_time=2.5e3;

[data, ic1,ic2] = par_col(ic1,ic2,N,max_time,fparam)