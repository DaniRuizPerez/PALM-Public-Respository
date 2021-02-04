classdef UnwrappedTSDBNSTS < StateTrackerSearch
% class for implementing network search algorithms that know which
% states need to be re-evalutated and which don't on Two-Stage Dynamic BNs
% Keeps lists of all possible edges and whether or not they are: legal,
% illegal, possible, and need to be reevaluated.
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.
    
    properties
        timeseries      % dynamic bayes nets
    end
    
    methods
        
        % constructor
        function dbnsts = UnwrappedTSDBNSTS(contData, discData, priorPrecision, ...
                phencol, backtracking, bf_thresh, nophenotype, checkRepeats, contNodeTypes, discNodeTypes)
            if (nargin < 5)
                backtracking = false;
            end
            if (nargin < 6)
                bf_thresh = 0;
            end
            if (nargin < 7)
                % change default here for dynamic bayes net:
                nophenotype = true;
            end
            if (nargin < 8)
                checkRepeats = true;
            end
            if (nargin < 9)
                contNodeTypes = [];
            end
            if (nargin < 10)
                discNodeTypes = [];
            end
            dbnsts@StateTrackerSearch(contData, discData, priorPrecision, ...
                phencol, backtracking, bf_thresh, nophenotype, checkRepeats, contNodeTypes, discNodeTypes);
            dbnsts.cycles = false;
            dbnsts.self = false;
            dbnsts.timeseries = true;
        end

        function obj = SetUnwrappedTimeSeries(obj)
            % call after dbnsts.Init()
        	% this sets up the search to be a dynamic bayesian network for
        	% time series data
            obj.timeseries = true;
            
            % unwrapped dynamic bayes nets care if there's cycles
            obj.cycles = false;
            
            % unwrapped dynamic bayes nets care if there's self-edges
            obj.self = false;
            
            % cut off edges from tn data to t0 data
            % continuous nodes are ordered before discrete nodes
            % and then doubled so we have :
            % [contdata_t0, contdata_tn, discdata_t0, discdata_tn]
            
            % close anything from contdata_t0 to contdata_t0:
            obj.closed(1:(obj.numContNode/2),1:(obj.numContNode/2)) = ones(obj.numContNode/2);
            % close anything from discdata_t0 to discdata_t0:
            obj.closed(obj.numContNode+1:obj.numContNode+(obj.numDiscNode/2), ...
                obj.numContNode+1:obj.numContNode+(obj.numDiscNode/2)) = ...
                ones(obj.numDiscNode/2);
            % close anything from discdata_t0 to contdata_t0:
            obj.closed(obj.numContNode+1:obj.numContNode+(obj.numDiscNode/2), ...
                1:(obj.numContNode/2)) = ...
                ones(obj.numDiscNode/2, obj.numContNode/2);            
            % close anything from discdata_tn to contdata_t0:
            obj.closed((obj.numContNode+1+(obj.numDiscNode/2)):end,1:(obj.numContNode/2)) = ...
                ones(obj.numDiscNode/2,obj.numContNode/2);
            % close anything from discdata_tn to discdata_t0
            obj.closed((obj.numContNode+1+(obj.numDiscNode/2)):end,(obj.numContNode+1):...
                (obj.numContNode+(obj.numDiscNode/2))) = ones(obj.numDiscNode/2);            
            % close anything from contdata_tn to contdata_t0
            obj.closed((obj.numContNode/2)+1:obj.numContNode,1:(obj.numContNode/2)) = ...
                ones(obj.numContNode/2);   
                        
            % update closed edges
            obj.never = obj.closed;
        end

        % check if edge is a self-edge 
        function selfEdge = isSelfEdge(obj, nodei, nodej)
            selfEdge = false;
            %continuous node
            if (nodei > 0 && nodei < obj.numContNode+1)
                if nodej == nodei + obj.numContNode/2
                    selfEdge = true;
                end
            %discrete node
            else
                if nodej == nodei + obj.numDiscNode/2
                    selfEdge = true;
                end
            end
        end
        
        function obj = SetInteractionMatrix(obj, interactionTypes, interactionMatrix)
            % override base class function 
            % intra-slice interactions between node types can be different in DBNs 
            % set up nodes restrictions (if any) from user-specified 
            % interaction matrix of node types 
            for i = 1:length(obj.bn)
                nodeIntrai = false;
                if (~ismember(obj.nodeDataTypes(i), interactionTypes))
                    continue
                end
                if ((i > obj.numContNode/2 && i < obj.numContNode+1) || (i > obj.numContNode+obj.numDiscNode/2))
                    nodeIntrai = true;
                end
                nodeTypei = find(interactionTypes == obj.nodeDataTypes(i));
                for j = 1:length(obj.bn)
                    nodeIntraj = false;
                    if (~ismember(obj.nodeDataTypes(j), interactionTypes) || i == j)
                        continue
                    end
                    if ((j > obj.numContNode/2 && j < obj.numContNode+1) || (j > obj.numContNode+obj.numDiscNode/2))
                        nodeIntraj = true;
                    end
                    nodeTypej = find(interactionTypes == obj.nodeDataTypes(j));
                    interactionType = interactionMatrix(nodeTypei, nodeTypej);
                    if (~interactionType)
                        obj.closed(i,j) = true;
                    elseif (~obj.isSelfEdge(i,j) && interactionType == 1)
                        obj.closed(i,j) = true;
                    elseif (nodeIntrai && nodeIntraj && interactionType == 2)
                        obj.closed(i,j) = true;
                    elseif (((nodeIntrai && ~nodeIntraj) || (~nodeIntrai && nodeIntraj)) && interactionType == 3)
                        obj.closed(i,j) = true;
                    end
                end
            end
            obj.never = obj.closed;
        end
        
        function [obj, numevals] = EvalSelfEdges(obj)
            % override base class function:
            % self-edges are indexed differently in DBNs with intra-slice edges
            numevals = 0;
            for i = 1:obj.numContNode/2
                % Consider candidate self-edge (i, i) for continuous node
                [obj, evalq] = obj.EvalOneEdge(i, i + obj.numContNode/2);
                if (evalq)
                    numevals = numevals + 1;
                end
            end
            for i = obj.numContNode+1:obj.numContNode+obj.numDiscNode/2
                % Consider candidate self-edge (i, i) for discrete node
                [obj, evalq] = obj.EvalOneEdge(i, i + obj.numDiscNode/2);
                if (evalq)
                    numevals = numevals + 1;
                end
            end
        end
        
    end
end
