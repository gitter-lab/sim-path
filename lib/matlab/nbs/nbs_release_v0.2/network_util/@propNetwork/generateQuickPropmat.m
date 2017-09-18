function obj = generateQuickPropmat(netobj,netpropFactor)
    
   obj = netobj;
   obj.propVal = (netpropFactor);
   
   %obj.adj_mat_norm_val = (1-netpropFactor)*inv(speye(size(netobj.adj_mat_norm)) - netpropFactor*netobj.adj_mat_norm); 
   
   obj.adj_mat_norm_val = (1-netpropFactor)*inv(eye(size(netobj.adj_mat_norm)) - netpropFactor*full(netobj.adj_mat_norm)); 
    
end