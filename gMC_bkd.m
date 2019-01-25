function sq=gMC_bkd(x)
% Input x array: (kf, km/kf, kt), where kf, kd, and kt represnt the bulge
% formaion, melting(dissolving), and traslocation rate.
%
% Matlab Parallel Computing Toolbox required.
% bulge formation rate (s-1)
kf=x(1);
% bulge melting rate (s-1)
km=x(2).*kf;
% bulge diffusion rate (s-1)
kt=x(3);
% Final State
finalStates=[6;9;11;12;13;14;15;16;17];

% Exp data
ExpData=[0.15649;0.21322;0.28986;0.30334;0.32609;0.3517;0.34632;0.36765;0.3632];
% sampling rate (s-1)
sampRate=1E6;
% Number of molecules (traces)
numMolec=5000;
% GPU operation
GPU=1;


NumFinalStates=size(finalStates,1);
fmtOutput=zeros(numMolec,NumFinalStates);
meanDwellTime=zeros(NumFinalStates,1);
if GPU==1
    unityArray=ones(numMolec*NumFinalStates,1,'gpuArray');
    parallel.gpu.rng();
else
    unityArray=ones(numMolec*NumFinalStates,1);
end

kfArray=kf*unityArray;
kmArray=km*unityArray;
kdArray=kt*unityArray;
for i=1:NumFinalStates
    finalStateArray((i-1)*numMolec+1:i*numMolec,1)=finalStates(i);
end

sampRateArray=sampRate*unityArray;

gOutput=arrayfun(@simuDwelltime,kfArray,kmArray,kdArray,finalStateArray,sampRateArray);
Output=gather(gOutput);
for i=1:NumFinalStates
    fmtOutput(:,i)=Output((i-1)*numMolec+1:i*numMolec);
    meanDwellTime(i,1)=mean(fmtOutput(:,i));
end

x
sq=sum((meanDwellTime-ExpData).^2)

end

function dwelltime=simuDwelltime(kf,km,kd,finalState,sampRate);
  t=0;
  state=0;
  while state~=finalState
      randNum=rand();
      t=t+1;
      if state==0
          if randNum<(kf/sampRate)
              state=1;
          else
              if randNum<(2*kf/sampRate)
                  state=-1;
              end
          end
      else
          if state==1
              if randNum<(km/sampRate)
                  state=0;
              else
                  if randNum<((km+kd)/sampRate)
                      state=2;
                  end
              end
          else
              if state==-1
                  if randNum<(km/sampRate)
                      state=0;
                  else
                      if randNum<((km+kd)/sampRate)
                          state=-2;
                      end
                  end
              else
                  if state==finalState*(-1)
                      if randNum<(kd/sampRate)
                          state=state+1;
                      end
                  else
                      if randNum<(kd/sampRate)
                          state=state+1;
                      else
                          if randNum<(2*kd/sampRate)
                              state=state-1;
                          end
                      end
                  end
              end
          end
      end   
  end
  dwelltime=t/sampRate;
end