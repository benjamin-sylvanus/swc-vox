classdef swc__vox
    %swc__vox 
    properties
        A;
        b;
        swc;
        R;
        batch;
    end

    methods
        function this = swc__vox(this)
%             addpath(genpath('/Users/benjaminsylvanus/Documents/GitHub/SparseMatrixGenerator'));
%             addpath(genpath('/Users/benjaminsylvanus/Documents/GitHub/plotcube'));
            this = this.init();
            this = this.batchInit()
            this = this.batchLoop()
        end

        function this = init(this)

            this.swc = swcreader();
            % Define Coordinates and Range
            X = this.swc.tree.X; Y = this.swc.tree.Y; Z = this.swc.tree.Z;
            dim = [X Y Z]; 
            
            % Threshold Radii
            [~,~,bin] = histcounts(this.swc.tree.Radii, 5);
            avg = mean(this.swc.tree.Radii(bin==1));
            this.swc.tree.Radii(this.swc.tree.Radii < avg)=avg;
            r = this.swc.tree.Radii;
            ranges = [min(dim - r); max(dim+r)]; % (min - ri) (max + ri)
            
            % Define Voxel Size
            VSIZE = avg/5;
            
            % Translate X,Y,Z
            tran = ((dim-ranges(1,:))./VSIZE)+1;
            
            % Update Bounds
            boundSize=ceil((ranges(2,:) - ranges(1,:))./VSIZE)+2;
            
            % Radii Scale
            this.R = this.swc.tree.Radii./(VSIZE);
            
            % Set Tree XYZ to Translated
            % min(this.swc.tree.Z-this.swc.tree.Radii);
            this.swc.tree.X = tran(:,1); this.swc.tree.Y = tran(:,2); this.swc.tree.Z = tran(:,3);
            this.swc.tree.Radii = this.R;
            % min(this.swc.tree.Z-this.swc.tree.Radii)
            
            this.b = this.calcBounds(this.swc,VSIZE);
            
            % Create Sparse
            this.A=ndSparse.build(boundSize);
        end


        function this = batchInit(this)
            swc = this.swc;
            tic;
            
            XD = swc.tree.X - swc.tree.X';
            YD = swc.tree.Y - swc.tree.Y';
            ZD = swc.tree.Z - swc.tree.Z';
            
            toc;
            
            tic;
            
            xds = XD.^2;
            yds = YD.^2;
            zds = ZD.^2;
            
            dists = sqrt(xds + yds + zds);
            
            toc;
            
            dist_id = dists < 50;
            full = [];
            this.batch = struct();
            c = 1;
            i = 1;
            while true
                temp = swc.tree(dist_id(:,i),:);
                a = ~ismember(temp.NodeId,full);
                full = [full; temp.NodeId(a)];
                rsd = swc.tree.NodeId(~ismember(swc.tree.NodeId,full));
                this.batch.Elements(c,1) = {temp.NodeId(a)};
                c=c+1;
                if numel(rsd>0)
                    i = rsd(1);
                else
                    break;
                end
            end
            
            % batch bound
            for i = 1:numel(this.batch.Elements)
                miX = 100000; miY = 100000; miZ = 100000;
                mxX = -1000; mxY = -1000; mxZ = -1000;
                ele = this.batch.Elements{i,1};
                comb = [miX miY miZ; mxX mxY mxZ];
                for j = 1:numel(this.batch.Elements{i,1})
                    Bound = this.b.a{ele(j),1};
                    Sx = Bound(1,1); Sy = Bound(1,2); Sz = Bound(1,3);
                    Nx = Bound(2,1); Ny = Bound(2,2); Nz = Bound(2,3);
                    small = find(Bound(1,:) < comb(1,:));
                    large = find(Bound(2,:) > comb(2,:));
            
                    comb(1,small) = Bound(1,small);
                    comb(2,large) = Bound(2,large);
                end
                this.batch.Bound(i,1) = {comb};
                this.batch.SIZE(i) = numel(ele);
            end
        end

        % Input [r1;r2] [x1 y1 z1; x2 y2 z2]
        % Output [min(ri-ci) max(ri-ci)]
        function b = calcBounds(this,swc,vs)    
            b = struct();
            t = swc.tree;
            mr = vs;
            for i = 1:numel(t.X)
                pair = swc.pairElement(i);
                b.a(i,:) =  this.pairBounds(t(i,:),pair,vs);         
            end
        
        end


        function b = pairBounds(this,p1, p2, m)
            % ri: pair radii
            % pxi: pair coords

            ri = [p1.Radii; p2.Radii]; 
            pxi = [p1{1,2:4}; p2{1,2:4}];
            b = {[round(min(pxi-ri)); ceil(max(pxi+ri))]};
            ba = [min(pi-(ri+m)), max(pi+(ri+m))];
        end


    % Replaced Calulation for r comparison
    function pos = swc2v(this,x0,y0,z0,x1,x2,y1,y2,z1,z2,r1,r2,Nx,Ny,Nz)
        t = ( (x0-x1)*(x2-x1) + (y0-y1)*(y2-y1) + (z0-z1)*(z2-z1) ) ./...
            ((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
        x = x1 + (x2-x1)*t;
        y = y1 + (y2-y1)*t;
        z = z1 + (z2-z1)*t;
    
        list1 = (x-x1).*(x-x2) + (y-y1).*(y-y2) + (z-z1).*(z-z2) <0;
        list2 = ~list1;
        
        dist2 = (x0-x).^2 + (y0-y).^2 + (z0-z).^2;
    
    
        %     r = r1 + sqrt((x-x1).^2 + (y-y1).^2 + (z-z1).^2) /...
        %         sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2) * (r2-r1);
    
        %     r = ( c + r2 ) / (sqrt ( 1 - ( |r1-r2 | / l ) )
    
        %     c = ( |r1 - r2| * l ) / L
    
            if r2 > r1 
                pos = swc2v(this,x0,y0,z0,x2,x1,y2,y1,z2,z1,r2,r1,Nx,Ny,Nz);
            else
        
        
            rd = abs( r1 - r2 );
        
            % distance from orthogonal vector to p2  
            l = sqrt( (x-x2).^2 + (y-y2).^2 + (z-z2).^2 );
        
            % distance from p1->p2     
            L = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);
        
            c = (rd * l) ./ sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);
            r = (c + r2) ./ sqrt(1 - ( (rd / L) .^2) );
        
        
            pos1 = dist2<=(r.^2); % smaller in one line and less than and equal
            pos2 = ( ( (x0-x1).^2 + (y0-y1).^2 + (z0-z1).^2 ) <= (r1^2) ) | ...
                   ( ( (x0-x2).^2 + (y0-y2).^2 + (z0-z2).^2 ) <= (r2^2) );
            pos = zeros(size(x0)); % use false
        
            
            pos(list1) = pos1(list1);
            pos(list2) = pos2(list2);
            end
    end



    function this = batchLoop(this)
        A=this.A;
        batch = this.batch;
        swc = this.swc;

        ecount = 0; prevElapse = 0; tic;

        for i = 1:numel(batch.Elements)

            Bound = this.batch.Bound{i,1};
    
            Sx = Bound(1,1) - 1; Sy = Bound(1,2) - 1; Sz = Bound(1,3) - 1;
    
            Nx = Bound(2,1) + 1; Ny = Bound(2,2) + 1; Nz = Bound(2,3) + 1;
    
            [y0,x0,z0]=meshgrid(Sy:Ny,Sx:Nx,Sz:Nz);
    
            pos_fill = ndSparse.build(range(Bound)+3);
    
            for j = 1:batch.SIZE(i)
    
                element = batch.Elements{i,1}(j);
    
                p1 = swc.tree{element,2:5}; 
    
                pair = swc.pairElement(element);
    
                p2 = pair{1,2:5};
    
                x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
    
                x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);
    
                pos = this.swc2v(x0,y0,z0,x1,x2,y1,y2,z1,z2,r1,r2,Nx,Ny,Nz);
    
                pos_fill = pos_fill | pos;
    
            end
    
            batch.TF(i) = {pos_fill};
    
            % To Plot:             
            %{
    
            hold on;
            is = isosurface(pos_fill,0);
            is.vertices = is.vertices + [Sy Sx Sz];
            p = patch('Faces',is.faces,'Vertices',is.vertices);
            p.FaceColor='green';
            p.EdgeColor='none';
            view(3);
            axis tight;
            daspect([1,1,.4]);
    
            %}
    
            ecount = ecount + this.batch.SIZE(i); 
    
            A(Sx:Nx,Sy:Ny,Sz:Nz) = pos_fill;
    
            telapse = toc;
    
            rate = this.batch.SIZE(i)/(telapse-prevElapse);
    
            prevElapse=telapse;
    
            s1 = sprintf(['I: %d \t Time: %f \t' ...
                ' Elements: %d \t Rate: %f \n'], ...
                        i,telapse, ecount,rate);
            disp(s1);
end

    end

    end

    
end

