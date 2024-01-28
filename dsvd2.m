function [normaerro2,cond2,maxc,atrab,matrizparaerro,matrizF,difM,vecsolsvd,novamatrizsistema1]=dsvd2(M,comeca,angulosft,anguloscol,anguloserro,listaint,rcolf,rcol,phicol,thetacol,thetacolf,phicolf,phierro,thetaerro,romega,vector,rerro,solucaoexacta)
quantostermos2= min((max(M,ceil(sqrt(size(angulosft,1)+1)-1))),15000);%listaint(end,1)%40; %%% número de termos na expansão do 'addition theorem'
quantostermos=(quantostermos2 +1)^(2);
matrizC=zeros(size(angulosft,1),quantostermos);
matrizF=zeros(quantostermos,size(anguloscol,1));
matrizparaerro=zeros(quantostermos,size(anguloserro,1));
%%%%%%%%%
for kk=1:quantostermos
matrizC(:,kk) = sphecomplex(listaint(kk,1), listaint(kk,2),phicolf,thetacolf)./(rcolf(:).^(listaint(kk,1)+1));
matrizF(kk,:) =((rcol(:).^(listaint(kk,1))).*conj(sphecomplex(listaint(kk,1),listaint(kk,2),phicol,thetacol)))./(romega*ones(size(rcol,1),1)).^(listaint(kk,1));%((raiocol(:).^(listaint(kk,1))).*///./(romega*ones(size(raiocol,1),1)).^(listaint(kk,1));
matrizparaerro(kk,:) = ((rerro(:).^(listaint(kk,1))).*conj(sphecomplex(listaint(kk,1),listaint(kk,2),phierro,thetaerro)))./(romega*ones(size(rerro,1),1)).^(listaint(kk,1));
end
MatD=zeros(size(matrizC,2));
for i=1:size(matrizC,2)
intcoluna=listaint(i,1);        
MatD(i,i)=(romega^(intcoluna))/((2*intcoluna+1));%(romega^(intcoluna))   4pi
end
difM=matrizC*MatD*matrizF;
atrab=matrizC*MatD;%*transpose(RR);%.*matrizFF);
[U,S,V]=svd(atrab,0);
V1=V(1:size(matrizC,1),:);
novamatrizsistema1=transpose(V1*matrizF);
matrizerro1=transpose(V1*matrizparaerro);  
ei=svd(novamatrizsistema1);
cond2=max(ei)/min(ei);
vecsolsvd=novamatrizsistema1 \ vector;
%%adaptar para obter erros da resolução do sistema.


maxc=max(vecsolsvd);
apr=matrizerro1*vecsolsvd;
erro=solucaoexacta-apr;
normaerro2=max(abs(erro))
end