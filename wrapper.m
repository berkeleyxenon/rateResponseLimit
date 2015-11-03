%

for source=1:4
typ=0;NRbandsim;plotResult;figure(501);save2pdf(['nu_case_' dis('%d',source) dis('_nco_%d',nco)]);
end