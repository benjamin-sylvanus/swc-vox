classdef swcreader
    %SWCREADER Summary of this class goes here
    %   Detailed explanation goes here

    properties
        filename;
        pathname;
        tree;
        pairs;
    end

    methods
        function this = swcreader(this)
            %SWCREADER Construct an instance of this class
            %   Detailed explanation goes here
            [this.filename, this.pathname] = uigetfile('*.swc','','');
            if isequal(this.filename,0)
                disp('User selected Cancel');
            else
                disp(['User selected ', fullfile(this.pathname,this.filename)]);
                this = this.read();
                this.pairElement(2);
            end

        end

        function this = read(this)
            fid = fopen(fullfile(this.pathname,this.filename));
            A = textscan (fid, '%s', 'delimiter', '\n');
            A = A{1};
            fclose (fid);
            swc = [];
            for counter  = 1 : length (A)
                if ~isempty (A{counter})  % allow empty lines in between
                    % allow comments: lines starting with #:
                    if ~strcmp (A{counter} (1), '#')
                        swc0   = sscanf (A{counter}, '%f')';
                        swc    = [swc; swc0];
                    end
                end
            end




            NodeID = swc(:,1); Coords = swc(:,3:5);
            Radii = swc(:,6); Parents = swc(:,7);

            this.tree = ...
                table(NodeID,...
                Coords(:,1),...
                Coords(:,2), ...
                Coords(:,3),...
                Radii, ...
                Parents, ...
                'VariableNames', ...
                {'NodeId', 'X','Y','Z', 'Radii', 'Parent'});
        end

        function pair = pairElement(this,id)
            %pairElement Selects the parent element of node defined by id.
            %   Returns parent element.
            %   Restricted to nodes(2:end); -- node id must be > 0.
            t= this.tree;
            parent = t.Parent(id);
            if parent~=-1
                pair = t(parent,:);
            else
                pair = t(1,:);
            end
            %             for i = 1:1000:numel(t.X)
            %                 parent = t.Parent(i);
            %                 if (parent~=-1)
            %                     disp('Node: ');
            %                     disp(t(i,:));
            %                     disp('Pair: ');
            %                     disp(t(parent,:));
            %                 end
            %             end
            %             disp(numel(this.tree.X))
        end
    end
end

