function [ LA,lambA,LB,lambB,U,e ] = itebdA( LA,lambA,LB,lambB,U,chimax,d,terror,tau )

        % applying the gate to the MPS wave function 
        psi=scon({lambB,LA,lambA,LB,lambB,U},{[-1 1],[1 5 2],[2 3],[ 3 6 4],[ 4 -4],[-2 -3 5 6 ]});
        %-log( trace(reshape(psi,chimax*d,chimax*d)*reshape(psi,chimax*d,chimax*d))  )/(0.001*2)

        gg=reshape(psi,chimax*d,chimax*d);

         % svding the time evolved state 
        [uu,dd,vv]=svd(reshape(psi,chimax*d,chimax*d),'econ');
        
       e=-log( scon({psi,psi},{[1 2 3 4] ,[1 2 3 4 ]}) )/tau/2;
       %gg-uu*dd*vv' 
       ttr=trace(dd*dd);
       
       % entanglement grows and bond dimesions should increase but I have to keep the computational resources fixed= truncation
 
       %======TRUNCATION=======
       uutrun=uu(1:chimax*d,1:chimax);
       ddtrun=dd(1:chimax,1:chimax);
       vv=vv';
       vvtrun=vv(1:chimax,1:chimax*d);
       error=1-trace(ddtrun*ddtrun)/ttr;
       ddtrun=ddtrun/sqrt(trace(ddtrun*ddtrun));
       lambA=ddtrun;
      
       if error>terror
           flag=1;
           tre=error;
           
           disp('large truncation')
           %i
           % disp(error)
           %break 
       end
       %====================== 
       
       %LA=reshape(uutrun,chimax,d,chimax); 
       %LB=reshape(vvtrun,chimax,d,chimax);
      
      
       [ invlb ] = invcare(lambB,chimax);
       % obtaining the time evolved wave function
       LB=reshape(reshape(vvtrun,chimax*d,chimax)*invlb,chimax,d,chimax);
       LA=reshape((invlb*reshape(uutrun,chimax,chimax*d)),chimax,d,chimax);
      

end

