function [id, unfID]= getIDpair(nid)

% get individual ID matrix
id = 1:nid;
ID = id'*ones(1,nid);
IDt = ID';
unfID = [ID(:) IDt(:)];