function [hat_x] = cs_result()
clc;clear
K=1;      
N=2000;    
M=2000;      
Phi=randn(M,N); 
a = Lorentzian()
x = [a(1:2000)]
s=Phi*x.';                                        
m=2*K;                                            
Psi=dwtmtx(2000,'coif5',3)./sqrt(N);                        
T=Phi*Psi';                                      
hat_y=zeros(1,N);                                                    
Aug_t=[];                                         
r_n=s;                                            
for times=1:1;                                    
    for col=1:N;                                  
        product(col)=abs(T(:,col)'*r_n);          
    end
    [val,pos]=max(product);                      
    Aug_t=[Aug_t,T(:,pos)];                       
    
    T(:,pos)=zeros(M,1);                          
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           
    r_n=s-Aug_t*aug_y;                            
    pos_array(times)=pos;                         
end
hat_y(pos_array)=aug_y; 
hat_x=abs(real(Psi'*hat_y.'));                         
figure(1);
plot(hat_x)
hold on;
plot(x,'r')
