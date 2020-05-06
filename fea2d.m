% This is the 2019 version of the FEA2d program for MAE-6010
function fea2d(varargin)
    fedata;
    ERROR_F = false;
    done = false;
    i=0;
    while(~done)
        i=i+1;
        % read in the input file name
        if(and(nargin>0,i<2))
            file1 = varargin{1};
        else
            file1 = input('Enter the input file name:   ','s'); %file1 = 'bridge.dat'; %
        end
        INFID = fopen(file1,'r');
        if(INFID < 0)
            disp('WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            disp('Invalid file name entered, try again');
        else
            done = true;
        end
    end
    if(nargin>1)
        file2 = varargin{2};
    else
        file2 = input('Enter the output file name:  ','s');  %file2 = 'bridge.out'; %
    end
    OUTFID = fopen(file2,'w');
    read_input;
    if(~ERROR_F) 
        assemble;  % solve for displacements and stresses
    end
    fclose('all');
    if(ERROR_F)
       disp('Errors occurred in the solution process.');
    else
       disp('Execution complete.');
    end
    i = input('Display Mesh? (1=yes): ');
    if(i==1)
        draw_mesh();
        j = input('Display Deformed Mesh? (1=yes): ');
        if(j==1)
            draw_deformed();
        end
    end
end 





function read_input()
    fedata;
    %  local variables
    a  = zeros(1,3);
    b  = zeros(1,3);
    c  = zeros(1,3);

    % initialize counters to zero and defaults for properties
    NUM_ELEM = 0;     % number of elements
    NUM_NODE = 0;     % number of nodes
    NUM_RBAR = 0;     % number of rigid bar elements read in
    NUM_SPC = 0;      % number of single point constraints
    NUM_LOADS = 0;    % number of loaded nodes
    NUM_MATERIAL = 0; % number of materials tables
    NUM_PHYS_PROP = 0;% number of physical property tables
    NUM_CS = 0;       % number of coordinate systems
    GX=0.;            % acceleration due to gravity in the x direction
    GY=0.;            % acceleration due to gravity in the y direction
    TEMP_INPUT = false; % flag to know if a temperature change was input
    title = blanks(80); % set model name to blank
    count_cards();

    frewind(INFID);        
    strng = blanks(80);
    done = false;
    % read the input file
    while(~done)
        strng = my_fgetl(INFID);
        if (strng == -1)  
            done = true; 
            break;
        end
        card = strtrim(strng(1:8)); % Renoves leading and trailing blanks
        ncc = length(card);
        card = [ card, blanks(8-ncc) ]; % Relace training blanks
        if( ncc == 0 )      % skip blank lines
            continue;
        elseif( strcmp(card(1:2),'ID') ) % ******** ID ****
            TITLE = strng(3:end);
            TITLE = strtrim(TITLE);
        elseif (strcmp(card(1:1),'$'))    % ******* $ ****
            continue;
        elseif (strcmp(card(1:4),'CBAR')) % ***************** CBAR ****
            NUM_ELEM = NUM_ELEM + 1;
            ne = str2num(strng( 9:16)); 
            pn = str2num(strng(17:24));
            n1 = str2num(strng(25:32));
            n2 = str2num(strng(33:40));
            ENODE(1,ne) = n1; 
            ENODE(2,ne) = n2;
            P_ID(ne) = pn;
            ELEM_TYPE(ne) = 2;
            % set z axis rotations free (constraints applied later)
            ID(6,n1) = 0;
            ID(6,n2) = 0;

            % **************************** CORD2R ****
        elseif (strcmp(card(1:6),'CORD2R') )
            NUM_CS = NUM_CS + 1;
            CS_ID(NUM_CS) = str2double(strng(9:16));
            rid = str2num(strng(17:24));
            a(1) = my_str2num(strng(25:32));
            a(2) = my_str2num(strng(33:40));
            a(3) = my_str2num(strng(41:48));
            b(1) = my_str2num(strng(49:56));
            b(2) = my_str2num(strng(57:64));
            b(3) = my_str2num(strng(65:72));
            if(rid ~= 0)
                disp('Error in a CORD2R input line.');
                disp('This program requires RID = 0');
                ERROR_F = true;
                return;
            end
            % we need to read another line
            if (check_plus(strng(73:80)))
                % this data is ignored
                strng = my_fgetl(INFID);
                if (strng == -1)
                    done = true;
                    disp('Missing 2nd line for CORD2R input')
                    ERROR_F = true;
                    break;
                end
                if (check_plus(strng(1:8)))
                    c(1) = my_str2num(strng( 9:16));
                    c(2) = my_str2num(strng(17:24));
                    c(3) = my_str2num(strng(25:32));
                else
                    disp('Error, missing line 2 for CORD2R card');
                    ERROR_F = true;
                    return;
                end
            else
                disp('Error, Missing line 2 for CORD2R card');
                ERROR_F = true;
                return;
            end
            % compute the transformation matrix and store in T_CS array
            make_cs_transform(a,b,c,NUM_CS);

            % ***************************** CQUAD4 ****
        elseif (strcmp(card(1:6),'CQUAD4'))
            NUM_ELEM = NUM_ELEM + 1;
            i      = str2double(strng( 9:16));
            P_ID(i)    = str2num(strng(17:24));
            ENODE(1,i) = str2num(strng(25:32));
            ENODE(2,i) = str2num(strng(33:40));
            ENODE(3,i) = str2num(strng(41:48));
            ENODE(4,i) = str2num(strng(49:56));
            ELEM_TYPE(i) = 4;

            % ***************************** CQUAD8 ****
        elseif (strcmp(card(1:6),'CQUAD8'))
            NUM_ELEM = NUM_ELEM + 1;
            i          = str2num(strng( 9:16));
            P_ID(i)    = str2num(strng(17:24));
            ENODE(1,i) = str2num(strng(25:32));
            ENODE(2,i) = str2num(strng(33:40));
            ENODE(3,i) = str2num(strng(41:48));
            ENODE(4,i) = str2num(strng(49:56));
            ENODE(5,i) = str2num(strng(57:64));
            ENODE(6,i) = str2num(strng(65:72));
            ELEM_TYPE(i) = 5;
            if (check_plus(strng(73:80)))  % we need to read another line
                strng = my_fgetl(INFID);
                if (strng == -1)
                    done = true;
                    break;
                end
                if(check_plus(strng(1:8)))
                    ENODE(7,i) = str2num(strng( 9:16));
                    ENODE(8,i) = str2num(strng(17:24));
                else
                    disp('Error, missing line 2 for CQUAD8 card');
                    ERROR_F = true;
                    return;
                end
            else
                disp('Error, missing line 2 for CQUAD8 card');
                ERROR_F = true;
                return;
            end

            % ******************************* CROD ****
        elseif (strcmp(card(1:4),'CROD'))
            NUM_ELEM = NUM_ELEM + 1;
            i          = str2num(strng( 9:16));
            P_ID(i)    = str2num(strng(17:24));
            ENODE(1,i) = str2num(strng(25:32));
            ENODE(2,i) = str2num(strng(33:40));
            ELEM_TYPE(i) = 1;

            % **************************** CTRIA3 ****
        elseif (strcmp(card(1:6),'CTRIA3'))
            NUM_ELEM = NUM_ELEM + 1;
            i          = str2num(strng( 9:16));
            P_ID(i)    = str2num(strng(17:24));
            ENODE(1,i) = str2num(strng(25:32));
            ENODE(2,i) = str2num(strng(33:40));
            ENODE(3,i) = str2num(strng(41:48));
            ELEM_TYPE(i) = 3;

            % ***************************** FORCE ****
        elseif (strcmp(card(1:5),'FORCE'))
            NUM_LOADS = NUM_LOADS + 1;
            LOADED_NODES(NUM_LOADS) = str2num(strng(17:24));
            f  = my_str2num(strng(33:40));
            fx = my_str2num(strng(41:48));
            fy = my_str2num(strng(49:56));
            F_X(NUM_LOADS)=f*fx;
            F_Y(NUM_LOADS)=f*fy;

            % ******************************* GRAV ****
        elseif (strcmp(card(1:4),'GRAV'))
            grav = my_str2num(strng(25:32));
            GX   = my_str2num(strng(33:40));
            GY   = my_str2num(strng(41:48));
            GX = grav*GX;
            GY = grav*GY;

            % ****************************** GRID ****
        elseif (strcmp(card(1:4),'GRID'))
            NUM_NODE = NUM_NODE + 1;
            i      = str2num(   strng( 9:16));
            X(i)   = my_str2num(strng(25:32));
            Y(i)   = my_str2num(strng(33:40));
            DCS(i) = str2num(   strng(49:56));
            cstr = strng(57:64);
            if(~isempty(strtrim(cstr)))
                i=0;
                store_spc(i, cstr, 0.);
                if(ERROR_F); return; end
            end

            % ****************************** MAT1 ****
        elseif (strcmp(card(1:4),'MAT1'))
            NUM_MATERIAL = NUM_MATERIAL + 1;
            i            = str2num(   strng( 9:16));
            E_MODULUS(i) = my_str2num(strng(17:24));
            G_MODULUS(i) = my_str2num(strng(25:32));
            PRATIO(i)    = my_str2num(strng(33:40));
            RHO(i)       = my_str2num(strng(41:48));
            ALPHA(i)     = my_str2num(strng(49:56));
            TREF(i)      = my_str2num(strng(57:64));
            if(PRATIO(i) == 0.)
                disp('A material with a blank or zero value for Poisson''s Ratio was found.  This is not allowed.');
                ERROR_F = true;
                return;
            end
            if(G_MODULUS(i) == 0.)
                G_MODULUS(i) = E_MODULUS(i)/( 2.*(1+PRATIO(i)) );
            end
            MAT_ID(NUM_MATERIAL) = i;
            % we need to read another line
            if(check_plus(strng(73:80)))
                % this data is ignored
                strng = my_fgetl(INFID);
                if (strng == -1)
                    done = true;
                    break;
                end
                if(check_plus(strng(1:8)))
                    continue;
                else
                    disp('Error, Missing line 2 for MAT1 card');
                    ERROR_F = true;
                    return;
                end
            end

            % ******************************* PBAR ****
        elseif (strcmp(card(1:4),'PBAR'))
            NUM_PHYS_PROP = NUM_PHYS_PROP + 1;
            i              = str2num(   strng( 9:16));
            M_ID(i)        = str2num(   strng(17:24));
            AREA(i)        = my_str2num(strng(25:32));
            MOM_INERTIA(i) = my_str2num(strng(33:40));
            PHYS_PROP_ID(NUM_PHYS_PROP) = i;
            % we need to read another line
            if(check_plus(strng(73:80)))
                strng = my_fgetl(INFID);
                if (strng == -1)
                    done = true;
                    break;
                end
                if(check_plus(strng(1:8)))
                    c1 = my_str2num(strng( 9:16));
                    c2 = my_str2num(strng(25:32));
                    c3 = my_str2num(strng(41:48));
                    c4 = my_str2num(strng(57:64));
                    H1(i) = max([c1,c2,c3,c4]);
                    H2(i) = min([c1,c2,c3,c4]);
                    % read another line
                    if(check_plus(strng(73:80)))
                        strng = my_fgetl(INFID);
                        if (strng == -1)
                            done = true;
                            break;
                        end
                        if(check_plus(strng(1:8)))
                            KSHEAR(i) = my_str2num(strng( 9:16));
                            if(KSHEAR(i) == 0.)
                                % a blank or zero value for shear ratio will
                                % be intrepreded to neglect the shear contribution
                                % a blank shear ratio
                                KSHEAR(i) = 1.0d38;
                            end
                        else
                            disp('Error, missing line 3 for PBAR card');
                            ERROR_F = true;
                            return;
                        end
                    else
                        KSHEAR(i) = 1.;
                    end
                else
                    disp('Error, missing line 2 for PBAR card');
                    ERROR_F = true;
                    return;
                end
            else
                disp('Error, missing line 2 for PBAR card');
                ERROR_F = true;
                return;
            end

            % ******************************* PROD ****
        elseif (strcmp(card(1:4),'PROD'))
            NUM_PHYS_PROP = NUM_PHYS_PROP + 1;
            i       = str2num(   strng( 9:16));
            M_ID(i) = str2num(   strng(17:24));
            AREA(i) = my_str2num(strng(25:32));
            PHYS_PROP_ID(NUM_PHYS_PROP) = i;

            % *************************** PSHELL ****
        elseif (strcmp(card(1:6),'PSHELL'))
            NUM_PHYS_PROP = NUM_PHYS_PROP + 1;
            i        = str2num(   strng( 9:16));
            M_ID(i)  = str2num(   strng(17:24));
            THICK(i) = my_str2num(strng(25:32));
            PHYS_PROP_ID(NUM_PHYS_PROP) = i;

            % ******************************* RBAR ****
        elseif (strcmp(card(1:4),'RBAR'))
            NUM_ELEM = NUM_ELEM + 1;
            NUM_RBAR = NUM_RBAR + 1;
            i          = str2num(strng( 9:16));
            ENODE(1,i) = str2num(strng(17:24));
            ENODE(2,i) = str2num(strng(25:32));
            str = strtrim(strng(33:40));
            ELEM_TYPE(i) = 6;
            RBAR_ID(NUM_RBAR) = i;
            % set z axis rotations free (constraints applied later)
            ID(6,ENODE(1,i)) = 0;
            ID(6,ENODE(2,i)) = 0;
            if ( ~strcmp(deblank(str),'123456'))
                disp(['Error: This program only supports RBAR elements', ...
                    'with all degrees of freedom slaved to the first node.']);
                ERROR_F = true;
                return;
            end

            % ******************************* SPC ****
        elseif (strcmp(card(1:4),'SPC '))
            i       = str2num(strng(17:24));
            cstr    = strng(25:32);
            value   = my_str2num(strng(33:40));
            store_spc(i, cstr, value);
            if(ERROR_F);  return; end

            % ****************************** SPC1 ****
        elseif (strcmp(card(1:4),'SPC1'))
            cstr = strng(17:24);
            g{1} = strng(25:32);
            g{2} = strng(33:40);
            g{3} = strng(41:48);
            g{4} = strng(49:56);
            g{5} = strng(57:64);
            g{6} = strng(65:72);
            g{2} = strtrim(g{2});
            % range of nodes specified
            if (strcmp(deblank(g{2}),'THRU'))
                i = str2num(g{1});
                j = str2num(g{3});
                for k = i: j
                    store_spc(k, cstr, 0.);
                    if(ERROR_F)
                        return;
                    end
                end
                k = j+1;
            else
                for i = 1: 6
                    k = str2num(g{i});
                    if(k ~= 0)
                        store_spc(k, cstr, 0.);
                        if(ERROR_F)
                            return;
                        end
                    else
                        break;
                    end
                end
            end

            % **************************** TEMPD ****
        elseif (strcmp(card(1:5),'TEMPD'))
            if(~TEMP_INPUT)
                TEMP_INPUT = true;
                TEMPD = my_str2num(strng(17:24));
            else
                disp('Error, multiple TEMPD cards read in.');
                ERROR_F = true;
                return;
            end

        else
            display(['Warning:  the following input card is ignored: ', strng]);
        end
    end

    % Set the constrained degrees of freedom in the ID matrix
    % to negative values
    for i = 1: NUM_SPC
        ID(SPC_DIR(i),SPC_NODE(i)) = -i;
    end

    % Add pseudo-nodes for the Lagrange multipliers for Rigid elements
    if(NUM_RBAR > 0)
        k = NUM_NODE;
        for i = 1: NUM_ELEM
            % it is an rbar element
            if(ELEM_TYPE(i) == 6)
                k = k + 1;
                % here we add a 3rd node to the RBAR element
                ENODE(3,i) = k;
            end
        end
    end

    return;
end % function read_input





function strng = my_fgetl(INFID)
    strng = fgetl(INFID);
    if (strng == -1) return;  end
    nc=length(strng);  % number of characters
    strng = [strng,blanks(80-nc)];  % Add trailing blanks
    return;
end




function val = my_str2num(str)
    val = 0.; 
    str = strtrim(str); % strip leading and trailing blanks
    L = length(str);
    ip = 0;
    ie = 0;
    for i = 1:L
        j = L + 1 - i;
        c = str(j:j);
        if(strcmp(c,'+')) 
            ip = j; 
        end
        if(strcmp(c,'-')) 
            ip = j; 
        end
        if(or(strcmp(c,'e'),strcmp(c,'E'))) 
            ie = j; 
        end
    end
    % Nastran numeric input does not include an E for the exponent to save
    % space.  The following code will insert one.
    if(and(ip > 1, ie == 0)) 
        str = [str(1:ip-1),'E',str(ip:L)];
    end
    if (L > 0)
        val = str2double(str);
    else
        val = 0.;
    end
end
    





function check_plusresult = check_plus(str)
    check_plusresult = false;
    str = strtrim(str); % remove leading blanks
    if(isempty(str))
        check_plusresult = false;
    elseif(strcmp(str(1:1),'+'))
        check_plusresult = true;
    end
    return;
end 






function store_spc(n, constr, value)
    fedata;

    % count restrained nodes and set restraints in the ID matrix
    for k = 1: 8
        dir = constr(k:k);
        if(strcmp(dir,'1'))
            NUM_SPC = NUM_SPC + 1; % count the # of constrainted dof
            SPC_NODE(NUM_SPC) = n; % store the node number
            SPC_DIR(NUM_SPC) = 1;  % store the direction
            CONSTR_DISP(NUM_SPC) = value; % store the displacement value
        elseif(strcmp(dir,'2'))
            NUM_SPC = NUM_SPC + 1;
            SPC_NODE(NUM_SPC) = n;
            SPC_DIR(NUM_SPC) = 2;
            CONSTR_DISP(NUM_SPC) = value;
        elseif(strcmp(dir,'3'))
            ID(3,n) = 1; % ignore this restraint in 2-dim. problems
        elseif(strcmp(deblank(dir),deblank('4')))
            ID(4,n) = 1; % ignore this restraint in 2-dim. problems
        elseif(strcmp(dir,'5'))
            ID(5,n) = 1; % ignore this restraint in 2-dim. problems
        elseif(strcmp(dir,'6'))
            NUM_SPC = NUM_SPC + 1;
            SPC_NODE(NUM_SPC) = n;
            SPC_DIR(NUM_SPC) = 6;
            CONSTR_DISP(NUM_SPC) = value;
        elseif(isempty(strtrim(dir)))
            continue;
        else
            str = sprintf('Unrecognized constraint for node =  %d',n);
            display(str);
            ERROR_F = true;
            return;
        end
    end
    return;
end % function store_spc






function assemble()
    fedata;
    % This program completes the following tasks
    %   1. assembles the stiffness matrix
    %   2. assembles the force vector
    %   3. sends the matrices for solution
    %   4. generates the displacement vectors
    %   5. calculates the element loads, strains and stresses
    % Method adapted from:
    %     Concepts and Applications of Finite Element Analysis
    %     3rd Edition
    %     Robert D Cook
    %     John Wiley, 1982, chapter 2

    % declare local variables
    kk  = int32(zeros(1,24));

    % Renumber the ID array
    % on input ID contains:
    %    zeros for d.o.f. in the problem
    %    ones for d.o.f. which are not in the solution process
    %    negative values for d.o.f. with single point constraints
    % on output ID contains:
    %    positive numbers for unknown d.o.f.
    %    zeros for d.o.f. which are not in the solution process
    %    negative values for d.o.f. with single point constraints
    % count the number of unknown displacements to be solved for
    N_UNKNOWN = 0;
    for n = 1: NUM_NODE
        for j = 1: 6
            if(ID(j,n) == 0 )
                N_UNKNOWN = N_UNKNOWN + 1;
                ID(j,n) = N_UNKNOWN;
            elseif(ID(j,n) == 1)
                ID(j,n) = 0;
            end
        end
    end

    % Now add the degrees of freedom associated with multipoint 
    % constraints from rigid elements. These will be the Lagrange Multipliers 
    % Essentially, each rigid element adds a new node to the ID array
    for n = NUM_NODE+1:NUM_NODE+NUM_RBAR
        N_UNKNOWN = N_UNKNOWN + 1;
        ID(1,n) = N_UNKNOWN;
        N_UNKNOWN = N_UNKNOWN + 1;
        ID(2,n) = N_UNKNOWN;
        ID(3,n) = 0;
        ID(4,n) = 0;
        ID(5,n) = 0;
        N_UNKNOWN = N_UNKNOWN + 1;
        ID(6,n) = N_UNKNOWN;
    end

    % Count the number of constrained d.o.f.
    ncons = 0;
    for n = 1: NUM_NODE
        for j = 1: 6
            if(ID(j,n) < 0)
                ncons = ncons + 1;
            end
        end
    end

    % Allocate arrays and set the contents S, R, and DISP arrays to zero
    R    = zeros(1,N_UNKNOWN+ncons);
    DISP = zeros(1,N_UNKNOWN+ncons);
    S    = zeros(N_UNKNOWN+ncons,N_UNKNOWN+ncons);

    % Now add the degrees of freedom associated with single point
    % constraints (enforced displacements) to the ID array.  These
    % d.o.f currently have negative values.
    % count the number of equations (unknown + known disp.)
    N_EQ = N_UNKNOWN;
    for n = 1: NUM_NODE
        for j = 1: 6
            if(ID(j,n) < 0)
                k = -ID(j,n);
                N_EQ = N_EQ + 1;
                ID(j,n) = N_EQ;
                % store the enforced displacement value
                DISP(N_EQ) = CONSTR_DISP(k);
            end
        end
    end

    % Assembly of stiffness matrix WITHOUT symmetric storage scheme.
    % This stiffness matrix contains both unknown and known d.o.f.
    mband = 1; % store the semibandwidth
    cnt   = 0; % counter for rbar elements
    n     = 0; % counter for element number
    for ii = 1: NUM_ELEM + NUM_RBAR
        n = ii;
        % Determine what kind of element it is and then compute the
        % element stiffness matrix and the locations where the terms
        % will be store in the structure stiffness matrix
        if(n <= NUM_ELEM)
            if(ELEM_TYPE(n) == 1)
                rod_elem(n);
                kk(1) = ID(1,ENODE(1,n));
                kk(2) = ID(2,ENODE(1,n));
                kk(3) = ID(1,ENODE(2,n));
                kk(4) = ID(2,ENODE(2,n));
                m = 4;
                nn = 2;
            elseif (ELEM_TYPE(n) == 2)
                bar_elem(n);
                kk(1) = ID(1,ENODE(1,n));
                kk(2) = ID(2,ENODE(1,n));
                kk(3) = ID(6,ENODE(1,n));
                kk(4) = ID(1,ENODE(2,n));
                kk(5) = ID(2,ENODE(2,n));
                kk(6) = ID(6,ENODE(2,n));
                m = 6;
                nn = 2;
            elseif (ELEM_TYPE(n) == 3)
                tria_elem(n);
                kk(1) = ID(1,ENODE(1,n));
                kk(2) = ID(2,ENODE(1,n));
                kk(3) = ID(1,ENODE(2,n));
                kk(4) = ID(2,ENODE(2,n));
                kk(5) = ID(1,ENODE(3,n));
                kk(6) = ID(2,ENODE(3,n));
                m = 6;
                nn = 3;
            elseif (ELEM_TYPE(n) == 4)
                if (USE_QM6)
                    qm6_elem(n);
                else
                    q4_elem(n);
                end;
                for i = 1:4;
                    j = 2*i;
                    kk(j-1) = ID(1,ENODE(i,n));
                    kk(j)   = ID(2,ENODE(i,n));
                end;
                m = 8;
                nn = 4;
            elseif (ELEM_TYPE(n) == 5)
                q8_elem(n);
                for i = 1:8
                    j = 2*i;
                    kk(j-1) = ID(1,ENODE(i,n));
                    kk(j)   = ID(2,ENODE(i,n));
                end;
                m = 16;
                nn = 8;
            elseif (ELEM_TYPE(n) == 6)  
                continue;  % process rbar elements later
            else
                str = sprintf('Element number %d  Element Type %d', n, ELEM_TYPE(n));
                disp('Error, an element currently not supported. ');
                disp(['    ',str]);
                ERROR_F = true;
                return;
            end
        else  % now process rbar elements
            cnt = cnt + 1;
            n = RBAR_ID(cnt);
            if (ELEM_TYPE(n) == 6)
                rbar_elem(n);
                kk(1) = ID(1,ENODE(1,n));
                kk(2) = ID(2,ENODE(1,n));
                kk(3) = ID(6,ENODE(1,n));
                kk(4) = ID(1,ENODE(2,n));
                kk(5) = ID(2,ENODE(2,n));
                kk(6) = ID(6,ENODE(2,n));
                kk(7) = ID(1,ENODE(3,n));
                kk(8) = ID(2,ENODE(3,n));
                kk(9) = ID(6,ENODE(3,n));
                m = 9;
                nn = 2;
            else
                display('Error, in adding rbar element stiffness terms.');
                str = sprintf('Element type = %d', ELEM_TYPE(n));
                ERROR_F = true;
                return;
            end
        end
        % Check that the ID array has returned a reasonable row or column location
        for i = 1: m;
            if(or(kk(i)<0 ,kk(i)>N_EQ))
                disp('Error in a reference to the ID array in function assemble');
                disp('You may have an element referencing an undefined node.');
                ERROR_F = true;
                if(ERROR_F)
                    break;
                end
            end
        end
        if(ERROR_F); return; end
        % If a node in the element stiffness matrix uses a displacement
        % coordinate system different that the base system (CS ID = 0),
        % apply a transformation of the stiffness terms associated with
        % that node.
        for i = 1: nn
            j = ENODE(i,n);
            if(DCS(j) ~= 0)  % new displacement C.S.
                if(ELEM_TYPE(n) ~= 6)  % do not transform rigid elements
                    trans_element(j, i, m, ELEM_TYPE(n));
                    if(ERROR_F)
                        return;
                    end
                end
            end
        end
        for i = 1:m
            k = kk(i);
            R(k) = R(k) + RE(i);
            for j = 1:m
                L = kk(j);
                S(k,L) = S(k,L) + SE(i,j);
                if(mband < L-k+1)
                    mband = L-k+1;
                end
            end
        end
    end

    % Assemble force vector
    for n = 1: NUM_LOADS
        i = LOADED_NODES(n);
        kx=ID(1,i);
        ky=ID(2,i);
        if(kx ~= 0)
            R(kx) = R(kx) + F_X(n);
        end
        if(ky ~= 0)
            R(ky) = R(ky) + F_Y(n);
        end
    end

    % If a node in the uses a displacement coordinate system different
    % that the base system (CS ID = 0), apply a transformation of the
    % force terms associated with that node.
    for n = 1: NUM_NODE
        % transformation needed.
        if(DCS(n) ~= 0)
            transform_force(n);
        end
    end

    % Find a small number that represents a numerical zero.
    % Check the diagonal terms of the stiffness matrix and find the 
    % largest term.  Using double precision math, a zero will be
    % assumed to be 12 digits smaller.
    eps = 0.;
    for i = 1:N_UNKNOWN
        if(abs(S(i,i)) > eps) 
            eps = abs(S(i,i));
        end
    end
    eps = eps*10.e-12;
    % Check for zero rotational dof stiffness associated with adding
    % rbar elements connected to plane stress elements.  We can
    % end up with a singular stiffness matrix if we don't add a
    % small stiffness for this situation.  The small stiffness we
    % will add is 10 times the ESPILON estimate of a zero.
    for ii = 1: NUM_RBAR
        n = RBAR_ID(ii);
        i = ID(6,ENODE(1,n));
        j = ID(6,ENODE(2,n));
        if (S(i,i) == 0.); S(i,i) = 10*eps; end;
        if (S(j,j) == 0.); S(j,j) = 10*eps; end;
    end;

    % Form a new load vector which includes the forces due to
    % known displacements using Eq. 2.7-3, page 40, in Cook, et. al.
    % 'Concepts and Applications of Finite Element Analysis, 4th Ed.'
    for i = 1: N_UNKNOWN
        for j = N_UNKNOWN + 1: N_EQ
            R(i) = R(i) - S(i,j)*DISP(j);
        end
    end

    % Solve [K]{D}={R}, displacements stored in the disp array
    DISP(1:N_UNKNOWN)=linsolve(S(1:N_UNKNOWN,1:N_UNKNOWN),R(1:N_UNKNOWN)');
    if(ERROR_F)
        return;
    end
    
    for i = N_UNKNOWN+1:N_EQ
        REACT_FORCE(i) = 0;
         for j = 1:N_EQ %N_UNKNOWN+1:N_EQ
            REACT_FORCE(i) = REACT_FORCE(i) + S(i,j)*DISP(j);
        end
        REACT_FORCE(i) = REACT_FORCE(i) - R(i);
    end
    
    
    % Generate the displacement vector for the whole truss
    % Displacements are stored in the disp array after solving
    for n = 1: NUM_NODE
        i = ID(1,n);
        j = ID(2,n);
        k = ID(6,n);
        if(i == 0 || j == 0)
            disp('Error!  Problem with ID array in resolving disp.');
            ERROR_F = true;
        end
        U(n) = DISP(i);
        V(n) = DISP(j);
        if(k == 0)
            % these d.o.f where not in the soltion.
            THETA_Z(n) = 0.;
        else
            THETA_Z(n) = DISP(k);
        end
    end
    % If a node uses a different displacement C.S., transform
    % displacements back to the global C.S.
    for n = 1: NUM_NODE
        if(DCS(n) ~= 0) % transformation needed.
            transform_disp(n);
        end
    end

    % Calculate loads and stresses in each element
    for n = 1: NUM_ELEM
        if(ELEM_TYPE(n) == 1)
            rod_stress(n);
        elseif (ELEM_TYPE(n) == 2)
            bar_stress(n);
        elseif (ELEM_TYPE(n) == 3)
            tria_stress(n);
        elseif (ELEM_TYPE(n) == 4)
            if (USE_QM6) 
                qm6_stress(n);
            else
                q4_stress(n);
            end
        elseif (ELEM_TYPE(n) == 5)
            q8_stress(n);
        end
    end

                       % Write out:
    disp_output();     % computed displacements
    rod_elem_out();    % Rod Element results
    bar_elem_out();  % Bar Element results
    tria_elem_out(); % Tria Element results
    q4_elem_out();   % Q4 Element results
    q8_elem_out();   % Q8 Element results
    reactf_output(); % Reaction Forces
    return;
end % function assemble







function rod_elem(n)
    fedata;

    idp = P_ID(n);
    idm = M_ID(idp);
    if(TEMP_INPUT)
        delta_t = TEMPD - TREF(idm);
    else
        delta_t = 0.;
    end
    n1 = ENODE(1,n);
    n2 = ENODE(2,n);
    L = sqrt((X(n2)-X(n1)).^2+(Y(n2)-Y(n1)).^2);
    sine =(Y(n2)-Y(n1))/L;
    cosine =(X(n2)-X(n1))/L;
    ael = AREA(idp)*E_MODULUS(idm)/L;
    c2  = cosine.^2*ael;
    cs  = cosine*sine*ael;
    s2  = sine.^2*ael;

    % Reallocate the element stiffness matrix and load vectors
    RE = zeros(4,1);
    SE = zeros(4,4);

    % Initialize the element stiffness matrix
    SE=[ c2   cs  -c2  -cs  ;  ...
         cs   s2  -cs  -s2  ;  ... 
        -c2  -cs   c2   cs  ;  ...
        -cs  -s2   cs   s2  ];

    % initialize the element force vector to the body forces
    mass = AREA(idp)*L*RHO(idm);
    RE(1) = GX*mass/2.;
    RE(2) = GY*mass/2.;
    RE(3) = RE(1);
    RE(4) = RE(2);

    % add forces due to temp. prestrains
    eaat = E_MODULUS(idm)*AREA(idp)*ALPHA(idm)*delta_t;
    RE(1) = RE(1) - cosine*eaat;
    RE(2) = RE(2) - sine*eaat;
    RE(3) = RE(3) + cosine*eaat;
    RE(4) = RE(4) + sine*eaat;

    return;
end % function rod_elem

function bar_elem(n)
    fedata;

    idp = P_ID(n);
    idm = M_ID(idp);
    if(TEMP_INPUT)
        delta_t = TEMPD - TREF(idm);
    else
        delta_t = 0.;
    end
    n1 = ENODE(1,n);
    n2 = ENODE(2,n);
    L = sqrt((X(n2)-X(n1)).^2+(Y(n2)-Y(n1)).^2);
    sine =(Y(n2)-Y(n1))/L;
    cosine =(X(n2)-X(n1))/L;
    area = AREA(idp);
    mom_iner = MOM_INERTIA(idp);
    G = G_MODULUS(idm);
    e = E_MODULUS(idm);
    k_y = 1/KSHEAR(idp);
    
    x = AREA(idp)*E_MODULUS(idm)/L;
    phi_y = 12*e*mom_iner*k_y/(area*G*L^2);
    y_1 = 12*e*mom_iner/((1+phi_y)*L^3);
    y_2 = 6*e*mom_iner/((1+phi_y)*L^2);
    y_3 = (4+phi_y)*e*mom_iner/((1+phi_y)*L);
    y_4 = (2-phi_y)*e*mom_iner/((1+phi_y)*L);
    
    f = x.*cosine.^2 + y_1*sine.^2;
    g = (x-y_1).*cosine.*sine;
    h = -y_2.*sine;
    q = y_2.*cosine;
    p = x*sine^2+y_1.*cosine^2;
 
    
%     c2  = cosine.^2*ael;
%     cs  = cosine*sine*ael;
%     s2  = sine.^2*ael;

    % Reallocate the element stiffness matrix and load vectors
    RE = zeros(6,1);
    SE = zeros(6,6);

    % Initialize the element stiffness matrix
    SE=[ f   g  h    -f   -g   h;  ...
         g   p  q    -g   -p   q;  ... 
         h   q  y_3  -h   -q   y_4;  ...
        -f  -g  -h    f    g   -h; ...
        -g  -p  -q    g    p   -q; ...
         h   q  y_4  -h    -q  y_3;];

    % initialize the element force vector to the body forces
    mass = AREA(idp)*L*RHO(idm);
    g_perp = -GX*sine+GY*cosine;
    RE(1) = GX*mass/2.;
    RE(2) = GY*mass/2.;
    RE(3) = mass*L*g_perp/12;
    RE(4) = RE(1);
    RE(5) = RE(2);
    RE(6) = -RE(3);
    % add forces due to temp. prestrains
    eaat = E_MODULUS(idm)*AREA(idp)*ALPHA(idm)*delta_t;
    RE(1) = RE(1) - cosine*eaat;
    RE(2) = RE(2) - sine*eaat;
    RE(4) = RE(4) + cosine*eaat;
    RE(5) = RE(5) + sine*eaat;

    return;
end % function bar_elem



% This function calculates the stiffness matrix for a constant
% strain triangle element (3 node).
function tria_elem(n)
    fedata;
    idp = P_ID(n);
    idm = M_ID(idp);
    if(TEMP_INPUT)
    delta_t = TEMPD - TREF(idm);
    else
    delta_t = 0.;
    end
    n1 = ENODE(1,n);
    n2 = ENODE(2,n);
    n3 = ENODE(3,n);
    x1 = X(n1);
    x2 = X(n2);
    x3 = X(n3);
    y1 = Y(n1);
    y2 = Y(n2);
    y3 = Y(n3);
    % calculate 2(area of triangle) and the area with absolute value
    a2 =(x2 - x1)*(y3 - y1) -(x3 - x1)*(y2 - y1);
    a = abs(a2)/2.;
    vol = a*THICK(idp);
    % Form [B]
    b=[(y2-y3)/a2     0.     (y3-y1)/a2     0.     (y1-y2)/a2    0.     ;  ...
    0.    (x3-x2)/a2      0.    (x1-x3)/a2     0.     (x2-x1)/a2 ;  ...
    (x3-x2)/a2 (y2-y3)/a2 (x1-x3)/a2 (y3-y1)/a2 (x2-x1)/a2 (y1-y2)/a2 ];
    % E matrix for plane stress, see Cook, page 21
    ee = E_MODULUS(idm)/(1. - PRATIO(idm)^2);
    en = ee*PRATIO(idm);
    G  = E_MODULUS(idm)/(2.*(1. + PRATIO(idm)));
    em = [ ee   en   0.  ;   ...
    en   ee   0.  ;   ...
    0.   0.   G  ] ;
    % Reallocate the element stiffness matrix and load vectors
    RE = zeros(6,1);
    SE = zeros(6,6);
    % form [B] transpose times [E]
    bteb = transpose(b)*em*b;
    % form [B]^t[E][B]tA
    SE(1:6,1:6) = bteb.*vol;
    % set element load vector to the body forces
    RE(1) = RHO(idm).*(GX*vol/3.);
    RE(2) = RHO(idm).*(GY*vol/3.);
    RE(3) = RE(1);
    RE(5) = RE(1);
    RE(4) = RE(2);
    RE(6) = RE(2);
    % add temperature prestrains to the element load vector
    eterm = E_MODULUS(idm)/(1. - PRATIO(idm));
    eat = eterm*ALPHA(idm)*delta_t*vol;
    RE(1) = RE(1) + eat*b(1,1);
    RE(2) = RE(2) + eat*b(2,2);
    RE(3) = RE(3) + eat*b(1,3);
    RE(4) = RE(4) + eat*b(2,4);
    RE(5) = RE(5) + eat*b(1,5);
    RE(6) = RE(6) + eat*b(2,6);
    return;
end % function tria_elem

function q4_elem(n)
   fedata;

    place  = [ 0.,    -0.577350269189626  -0.774596669241483; ...
               0.,     0.577350269189626,  0.; ...
               0.,     0.,                 0.774596669241483 ];
    wgt =  [ 2.,      1.,     0.555555555555556; ...
             0.,      1.,     0.888888888888889; ... 
             0.,      0.,     0.555555555555556 ];

    % define local variables
    % em(i,j) = modulus matrix
    % thc(i) = nodal thicknesses
    % b(i,j) = strain displacement matrix
    % bte(i,j) = transpose of b(i,j) times em(i,j)
    % place(i,j) = sampling points for gaussian quadrature
    %              j = order of quadrature to be used (nguass=1,2,3)
    %              i = 1, 2, or 3 for the sampling points
    % xi = a particular gaussian quadrature sampling point
    % eta = a particular gaussian quadrature sampling point
    % wgt(i,j) = wieghts for gaussian quadrature
    %              j = order of quadrature to be used (nguass=1,2,3)
    %              i = 1, 2, or 3 for the weight
    % n(i) = shape function ni, i=1, 2, 3, or 4
    %        n(i) is computed by subroutine shape at specified values
    %        of xi and eta
    % exi(i) = initial strain in the x direction at node i, i=1,2,3,4
    % eyi(i) = initial strain in the y direction at node i, i=1,2,3,4
    % esi(i) = initial shear strain at node i, i=1,2,3,4
    % eth(i) = interpolated strains at a gauss point
    % thc(i) = thickness of element at node i, i=1,2,3,4
    % thk = thickness of element at specified values of xi and eta

    % number of nodes in the element
    nn = 4;
    % number of degrees of freedom in the element
    nd = nn*2;
    idp = P_ID(n);
    idm = M_ID(idp);
    if(TEMP_INPUT)
        delta_t = TEMPD - TREF(idm);
    else
        delta_t = 0.;
    end

    % store node numbers and their x and y coordinates
    for i = 1: 4
        xl(i)  = X(ENODE(i,n));
        yl(i)  = Y(ENODE(i,n));
    end

    % initialize the nodal thickness and initial strains
    for i=1:4
        thc(i) = THICK(idp);
        exi(i) = ALPHA(idm)*delta_t;
        eyi(i) = ALPHA(idm)*delta_t;
        esi(i) = 0.;
    end

    % initialize the modulus matrix
    % E matrix for plane stress, see Cook, page 21
    ee = E_MODULUS(idm)/(1. - PRATIO(idm)^2);
    en = ee*PRATIO(idm);
    G  = E_MODULUS(idm)/(2.*(1. + PRATIO(idm)));
    em = [ ee   en   0.  ;   ...
           en   ee   0.  ;   ...
           0.   0.   G  ] ;

    % Reallocate the element stiffness matrix and load vectors
    RE = zeros(8,1);
    SE = zeros(8,8);

    % start guass quadrature loop.  use ngauss by ngauss rule.
    ngauss = 2;
    for na = 1:ngauss
        xi = place(na,ngauss);
        for nb = 1:ngauss
            eta = place(nb,ngauss);
            [eth,thk,detjac,ns,b]=q4_shape(xi,eta,xl,yl,exi,eyi,esi,thc);
            dv = abs(wgt(na,ngauss)*wgt(nb,ngauss)*thk*detjac);
            % store b-transpose times e in 8 by 3 work-space array bte.
            for j = 1: nn
                L = 2*j;
                k = L-1;
                % do only multiplications that give a nonzero product.
                for m = 1:3
                    bte(k,m) = b(1,k)*em(1,m) + b(3,k)*em(3,m);
                    bte(L,m) = b(2,L)*em(2,m) + b(3,L)*em(3,m);
                end
                % add contribution of body forces to element nodal load array
                RE(k) = RE(k) + ns(j)*RHO(idm)*GX*dv;
                RE(L) = RE(L) + ns(j)*RHO(idm)*GY*dv;
            end
            % loop on rows of arrays se and re
            for nrow = 1: nd
                % add contribution of initial strains to load array re.
                for j = 1: 3
                    RE(nrow) = RE(nrow) + bte(nrow,j)*eth(j)*dv;
                end
                % loop to add contribution to element stiffness array se.
                for ncol = nrow: nd
                    dum=0.;
                    %  loop for product (b)t*(e)*(b).  zeros in b not skipped.
                    for j = 1: 3
                        dum = dum + bte(nrow,j)*b(j,ncol);
                    end
                    SE(nrow,ncol) = SE(nrow,ncol) + dum*dv;
                end
            end
        end
    end

    % fill in lower triangle of element stiffness matrix by symmetry.
    for i = 1: nd-1
        for j = i: nd
            SE(j,i)= SE(i,j);
        end
    end

    return;
end % function q4_elem

function [eth,thk,detjac,ns,b]=q4_shape(xi,eta,xl,yl,exi,eyi,esi,thc)
    % This program was adapted from:
    % 'Concepts and Application of Finite Element Analysis'
    % Robert D Cook, John Wiley, 1989, page 174, figure 6.5-1
    % define local variables
    % Nxi(i) = derivative of n(i) with respect to xi
    % Net(i) = derivative of n(i) with respect to eta
    % jac(i,j) = jacobian matrix
    % Shape Functions
    ns(1) =(1. - xi)*(1. - eta)/4.;
    ns(2) =(1. + xi)*(1. - eta)/4.;
    ns(3) =(1. + xi)*(1. + eta)/4.;
    ns(4) =(1. - xi)*(1. + eta)/4.;
    % Derivative of shape functions with respect to xi
    nxi(1) = -(1. - eta)/4.;
    nxi(2) =  (1. - eta)/4.;
    nxi(3) =  (1. + eta)/4.;
    nxi(4) = -(1. + eta)/4.;
    % Derivative of shape functions with respect to eta
    net(1) = -(1. - xi)/4.;
    net(2) = -(1. + xi)/4.;
    net(3) =  (1. + xi)/4.;
    net(4) =  (1. - xi)/4.;
    % clear arrays jac, eth and b
    jac = zeros(2,2);
    % Find jacobian jac and its determinant.  Replace jac by its inverse
    for i=1:4
    jac(1,1) = jac(1,1) + nxi(i)*xl(i);
    jac(1,2) = jac(1,2) + nxi(i)*yl(i);
    jac(2,1) = jac(2,1) + net(i)*xl(i);
    jac(2,2) = jac(2,2) + net(i)*yl(i);
    end
    detjac = jac(1,1)*jac(2,2) - jac(2,1)*jac(1,2);
    tmp = jac(1,1)/detjac;
    jac(1,1) =  jac(2,2)/detjac;
    jac(1,2) = -jac(1,2)/detjac;
    jac(2,1) = -jac(2,1)/detjac;
    jac(2,2) =  tmp;
    % form strain-displacement matrix b (zero entries are already set)
    b = zeros(3,8);
    for j = 1: 4
    L = 2*j;
    k = L - 1;
    b(1,k) = jac(1,1)*nxi(j) + jac(1,2)*net(j);
    b(2,L) = jac(2,1)*nxi(j) + jac(2,2)*net(j);
    b(3,k) = b(2,L);
    b(3,L) = b(1,k);
    end
    % interpolate initial strains and thickness from corner values
    thk = 0.;
    eth = zeros(1,3);
    for i=1:4
    eth(1) = eth(1) + ns(i)*exi(i);
    eth(2) = eth(2) + ns(i)*eyi(i);
    eth(3) = eth(3) + ns(i)*esi(i);
    thk = thk + ns(i)*thc(i);
    end
    return;
end % function q4_shape

function qm6_elem(n)
    fedata;

    place  = [ 0.,    -0.577350269189626  -0.774596669241483; ...
               0.,     0.577350269189626,  0.; ...
               0.,     0.,                 0.774596669241483 ];
    wgt =  [ 2.,      1.,     0.555555555555556; ...
             0.,      1.,     0.888888888888889; ... 
             0.,      0.,     0.555555555555556 ];

    % define local variables
    % em(i,j) = modulus matrix
    % thc(i) = nodal thicknesses
    % b(i,j) = strain displacement matrix
    %           note: the last column in b(i,j) holds the interpolated
    %                values of initial strains in the x and y directions
    %                and the shear strain
    % bte(i,j) = transpose of b(i,j) times em(i,j)
    % place(i,j) = sampling points for gaussian quadrature
    %              j = order of quadrature to be used (nguass=1,2,3)
    %              i = 1, 2, or 3 for the sampling points
    % xi = a particular gaussian quadrature sampling point
    % eta = a particular gaussian quadrature sampling point
    % wgt(i,j) = wieghts for gaussian quadrature
    %              j = order of quadrature to be used (nguass=1,2,3)
    %              i = 1, 2, or 3 for the weight
    % ns(i) = shape function ni, i=1, 2, 3, or 4
    %        ns(i) is computed by subroutine shape at specified values
    %        of xi and eta
    % exi(i) = initial strain in the x direction at node i, i=1,2,3,4
    % eyi(i) = initial strain in the y direction at node i, i=1,2,3,4
    % esi(i) = initial shear strain at node i, i=1,2,3,4
    % eth(i) = interpolated strains at a gauss point
    % thc(i) = thickness of element at node i, i=1,2,3,4
    % thk = thickness of element at specified values of xi and eta

    % number of nodes in the element (counting nodeless d.o.f.)
    nn  = 6;
    % number of degrees of freedom in the element
    nd  = nn*2;
    idp = P_ID(n);
    idm = M_ID(idp);
    if(TEMP_INPUT)
        delta_t = TEMPD - TREF(idm);
    else
        delta_t = 0.;
    end

    % store node numbers and their x and y coordinates
    for i = 1: 4
        xl(i)  = X(ENODE(i,n));
        yl(i)  = Y(ENODE(i,n));
    end

    % initialize the nodal thickness and initial strains
    for i=1:4
        thc(i) = THICK(idp);
        exi(i) = ALPHA(idm)*delta_t;
        eyi(i) = ALPHA(idm)*delta_t;
        esi(i) = 0.;
    end

    % initialize the modulus matrix
    % E matrix for plane stress, see Cook, page 21
    ee = E_MODULUS(idm)/(1. - PRATIO(idm)^2);
    en = ee*PRATIO(idm);
    G  = E_MODULUS(idm)/(2.*(1. + PRATIO(idm)));
    em = [ ee   en   0.  ;   ...
           en   ee   0.  ;   ...
           0.   0.   G  ] ;

    % clear element load vector and the stiffness mat.
    RE = zeros(12,1);
    SE = zeros(12,12);

    % start guass quadrature loop.  use ngauss by ngauss rule.
    ngauss = 2;
    for na = 1:ngauss
        xi = place(na,ngauss);
        for nb = 1:ngauss
            eta = place(nb,ngauss);
            [eth,thk,detjac,ns,b,djz]=qm6_shape(xi,eta,xl,yl,exi,eyi,esi,thc);
            dv  = abs(wgt(na,ngauss)*wgt(nb,ngauss)*thk*detjac);
            dvz = abs(wgt(na,ngauss)*wgt(nb,ngauss)*thk*djz);
            % Store [B]-transpose times [E] in 12 by 3 work-space array bte.
            for j = 1: nn
                L = 2*j;
                k = L-1;
                % Do only multiplications that give a nonzero product.
                for m = 1:3
                    bte(k,m) = b(1,k)*em(1,m) + b(3,k)*em(3,m);
                    bte(L,m) = b(2,L)*em(2,m) + b(3,L)*em(3,m);
                end
                % Add contribution of body forces to element nodal load array
                % ignore body force of nodeless dof
                if(j < 5)
                    RE(k) = RE(k) + ns(j)*RHO(idm)*GX*dv;
                    RE(L) = RE(L) + ns(j)*RHO(idm)*GY*dv;
                end
            end
            % Loop on rows of arrays se and re
            for nrow = 1: nd
                % Add contribution of initial strains to load array re.
                for j = 1: 3
                    % ignore thermal force for nodeless dof
                    if(nrow < 9)
                        RE(nrow) = RE(nrow) + bte(nrow,j)*eth(j)*dv;
                    end
                end
                % Loop to add contribution to element stiffness array se.
                for ncol = nrow: nd
                    dum = 0.;
                    %  Loop for product (b)t*(e)*(b).  Zeros in [B] not skipped.
                    for j = 1: 3
                        dum = dum + bte(nrow,j)*b(j,ncol);
                    end
                    if(or(nrow > 8, ncol > 8))
                        SE(nrow,ncol) = SE(nrow,ncol) + dum*dvz;
                    else
                        SE(nrow,ncol) = SE(nrow,ncol) + dum*dv;
                    end
                end
            end
        end
    end

    % Fill in lower triangle of element stiffness matrix by symmetry.
    for i = 1: nd-1
        for j = i: nd
            SE(j,i)= SE(i,j);
        end
    end

    % Condense out the nodeless d.o.f
    condense(4,12);
end % function qm6_elem

function [eth,thk,detjac,ns,b,detjacz] = qm6_shape(xi,eta,xl,yl,exi,eyi,esi,thc)

    % This program was adapted from:
    % 'Concepts and Application of Finite Element Analysis'
    % Robert D Cook, John Wiley, 1989, page 174, figure 6.5-1

    % define local variables
    % Nxi(i) = derivative of n(i) with respect to xi
    % Net(i) = derivative of n(i) with respect to eta
    % jac(i,j) = jacobian matrix
    % xii(i) = constants used in calculating the shape function
    % eti(i) = constants used in calculating the shape function

    % Shape Functions
    ns  = zeros(1,6);
    ns(1) =(1. - xi)*(1. - eta)/4.;
    ns(2) =(1. + xi)*(1. - eta)/4.;
    ns(3) =(1. + xi)*(1. + eta)/4.;
    ns(4) =(1. - xi)*(1. + eta)/4.;
    ns(5) = 1. - xi^2;
    ns(6) = 1. - eta^2;
    % Derivative of shape functions with respect to xi
    nxi(1) = -(1. - eta)/4.;
    nxi(2) =  (1. - eta)/4.;
    nxi(3) =  (1. + eta)/4.;
    nxi(4) = -(1. + eta)/4.;
    nxi(5) = -2*xi;
    nxi(6) = 0.;
    % Derivative of shape functions with respect to eta
    net(1) = -(1. - xi)/4.;
    net(2) = -(1. + xi)/4.;
    net(3) =  (1. + xi)/4.;
    net(4) =  (1. - xi)/4.;
    net(5) = 0.;
    net(6) = -2*eta;

    % Clear arrays jac
    jac = zeros(2,2);
    % Find jacobian jac and its determinant.  Replace jac by its inverse.
    for i = 1: 4
        jac(1,1) = jac(1,1) + nxi(i)*xl(i);
        jac(1,2) = jac(1,2) + nxi(i)*yl(i);
        jac(2,1) = jac(2,1) + net(i)*xl(i);
        jac(2,2) = jac(2,2) + net(i)*yl(i);
    end
    detjac = jac(1,1)*jac(2,2) - jac(2,1)*jac(1,2);
    tmp = jac(1,1)/detjac;
    jac(1,1) =  jac(2,2)/detjac;
    jac(1,2) = -jac(1,2)/detjac;
    jac(2,1) = -jac(2,1)/detjac;
    jac(2,2) =  tmp;

    % Find the jacobian and its determinant with xi = eta = 0 for use in
    % forming the terms last 4 columns to make a QM6 element
    jacz(1,1) =(-xl(1) + xl(2) + xl(3) - xl(4))/4.;
    jacz(1,2) =(-yl(1) + yl(2) + yl(3) - yl(4))/4.;
    jacz(2,1) =(-xl(1) - xl(2) + xl(3) + xl(4))/4.;
    jacz(2,2) =(-yl(1) - yl(2) + yl(3) + yl(4))/4.;
    detjacz = jacz(1,1)*jacz(2,2) - jacz(2,1)*jacz(1,2);
    tmp = jacz(1,1)/detjacz;
    jacz(1,1) =  jacz(2,2)/detjacz;
    jacz(1,2) = -jacz(1,2)/detjacz;
    jacz(2,1) = -jacz(2,1)/detjacz;
    jacz(2,2) =  tmp;

    % Form strain-displacement matrix [B](zero entries are already set).
    % This will form colums 1 through 8
    b   = zeros(3,10);
    for j = 1: 4
        L = 2*j;
        k = L - 1;
        b(1,k) = jac(1,1)*nxi(j) + jac(1,2)*net(j);
        b(2,L) = jac(2,1)*nxi(j) + jac(2,2)*net(j);
        b(3,k) = b(2,L);
        b(3,L) = b(1,k);
    end
    % Now form the last four columns of [B]
    for j = 5: 6
        L = 2*j;
        k = L - 1;
        b(1,k) = jacz(1,1)*nxi(j) + jacz(1,2)*net(j);
        b(2,L) = jacz(2,1)*nxi(j) + jacz(2,2)*net(j);
        b(3,k) = b(2,L);
        b(3,L) = b(1,k);
    end

    % Interpolate initial strains and thickness from corner values.
    % Note that a bilinear interpolation is used.
    thk = 0.;
    eth = zeros(1,3);
    for i=1: 4
        eth(1) = eth(1) + ns(i)*exi(i);
        eth(2) = eth(2) + ns(i)*eyi(i);
        eth(3) = eth(3) + ns(i)*esi(i);
        thk = thk + ns(i)*thc(i);
    end
    return;
end % function qm6_shape


function condense(num_condense,nsize)
    fedata;

    % The algorithm here is from Cook, 'Concepts and Applications of
    % Finite Element Analysis', 3 Ed. page 231
    % do condensation over the lower triangle
    for k = 1: num_condense
        n = nsize - k;
        j = n + 1;
        for i = 1: n
            tmp = SE(j,i)/SE(j,j);
            for m = 1: i
                SE(i,m) = SE(i,m) - SE(j,m)*tmp;
            end; m = i+1;
            RE(i) = RE(i) - RE(j)*tmp;
        end
    end
    % now fill in the upper triangle of se by symmetry
    for j = 1: n
        for i = 1: j
            SE(i,j) = SE(j,i);
        end
    end
end % function condense

function q8_elem(n)
    % This function computes the element stiffness matrix for a quadratic
    % quadrilateral plane element assuming plane stress loading.  It forms
    % a Lagrange type element (9 node) and condenses the center node.
    fedata;

    place  = [ 0.,    -0.577350269189626  -0.774596669241483; ...
               0.,     0.577350269189626,  0.; ...
               0.,     0.,                 0.774596669241483 ];
    wgt =  [ 2.,      1.,     0.555555555555556; ...
             0.,      1.,     0.888888888888889; ... 
             0.,      0.,     0.555555555555556 ];

    % define local variables
    % em(i,j) = modulus matrix
    % thc(i) = nodal thicknesses
    % b(i,j) = strain displacement matrix
    %           note: the last column in b(i,j) holds the interpolated
    %                values of initial strains in the x and y directions
    %                and the shear strain
    % bte(i,j) = transpose of b(i,j) times em(i,j)
    % place(i,j) = sampling points for gaussian quadrature
    %              j = order of quadrature to be used (nguass=1,2,3)
    %              i = 1, 2, or 3 for the sampling points
    % xi = a particular gaussian quadrature sampling point
    % eta = a particular gaussian quadrature sampling point
    % wgt(i,j) = wieghts for gaussian quadrature
    %              j = order of quadrature to be used (nguass=1,2,3)
    %              i = 1, 2, or 3 for the weight
    % n(i) = shape function ni, i=1, 2, 3, or 4
    %        n(i) is computed by subroutine shape at specified values
    %        of xi and eta
    % exi(i) = initial strain in the x direction at node i, i=1,2,3,4, ... 8
    % eyi(i) = initial strain in the y direction at node i, i=1,2,3,4, ... 8
    % esi(i) = initial shear strain at node i, i=1,2,3,4, ... 8
    % eth(i) = interpolated strains at a gauss point
    % thc(i) = thickness of element at node i, i=1,2,3,4, ... 8
    % thk = thickness of element at specified values of xi and eta

    % number of nodes in the element
    nn = 9;
    % number of degrees of freedom in the element
    nd = nn*2;
    idp = P_ID(n);
    idm = M_ID(idp);
    if(TEMP_INPUT)
        delta_t = TEMPD - TREF(idm);
    else
        delta_t = 0.;
    end;

    % get the attached node numbers and coordinates
    for i = 1: nn-1
        xl(i) = X(ENODE(i,n));
        yl(i) = Y(ENODE(i,n));
    end;

    % initialize the nodal thickness and initial strains
    for i=1:nn
        thc(i) = THICK(idp);
        exi(i) = ALPHA(idm)*delta_t;
        eyi(i) = ALPHA(idm)*delta_t;
        esi(i) = 0.;
    end;

    % initialize the modulus matrix
    % E matrix for plane stress, see Cook, page 21
    ee = E_MODULUS(idm)/(1. - PRATIO(idm)^2);
    en = ee*PRATIO(idm);
    G  = E_MODULUS(idm)/(2.*(1. + PRATIO(idm)));
    em = [ ee   en   0.  ;   ...
           en   ee   0.  ;   ...
           0.   0.   G  ] ;

    % clear element load vector and upper triangle of el. stiffness mat.
    RE = zeros(18,1);
    SE = zeros(18,18);

    % start guass quadrature loop.  use ng by ng gauss quadrature rule.
    ng = 3;
    for na = 1:ng
        xi = place(na,ng);
        for nb = 1:ng
            eta = place(nb,ng);
            [eth,thk,detjac,ns,b]=q8_shape(xi,eta,xl,yl,exi,eyi,esi,thc);
            dv = abs(wgt(na,ng)*wgt(nb,ng)*thk*detjac);
            % store b-transpose times e in 18 by 3 work-space array bte.
            for j = 1: nn
                L = 2*j;
                k = L-1;
                % do only multiplications that give a nonzero product.
                for m = 1:3
                    bte(k,m) = b(1,k)*em(1,m) + b(3,k)*em(3,m);
                    bte(L,m) = b(2,L)*em(2,m) + b(3,L)*em(3,m);
                end;
                % add contribution of body forces to element nodal load array
                RE(k) = RE(k) + ns(j)*RHO(idm)*GX*dv;
                RE(L) = RE(L) + ns(j)*RHO(idm)*GY*dv;
            end;
            % loop on rows of arrays SE and re
            for nrow = 1: nd
                % add contribution of initial strains to load array re.
                for j = 1: 3
                    RE(nrow) = RE(nrow) + bte(nrow,j)*eth(j)*dv;
                end;
                % loop to add contribution to element stiffness array SE.
                for ncol = nrow: nd
                    dum = 0.;
                    %  loop for product (b)t*(e)*(b).  zeros in b not skipped.
                    for j = 1: 3
                        dum = dum + bte(nrow,j)*b(j,ncol);
                    end; 
                    SE(nrow,ncol) = SE(nrow,ncol) + dum*dv;
                end;
            end;
        end;

    end;

    % fill in lower triangle of element stiffness matrix by symmetry.
    for i = 1: nd - 1
        for j = i: nd
            SE(j,i)= SE(i,j);
        end;
    end; 

    % Condense out the center d.o.f
    condense(2,18);

    return;
end % function q8_elem






function [eth,thk,detjac,ns,b]=q8_shape(xi,eta,xl,yl,exi,eyi,esi,thc)
    % Define local variables
    % Nxi(i) = derivative of n(i) with respect to xi
    % Net(i) = derivative of n(i) with respect to eta
    % jac(i,j) = jacobian matrix

    % Using serendipity element formulation, find x and y coordiantes
    % for a center node
    xl(9) = -(xl(1)+xl(2)+xl(3)+xl(4))/4. +(xl(5)+xl(6)+xl(7)+xl(8))/2.;
    yl(9) = -(yl(1)+yl(2)+yl(3)+yl(4))/4. +(yl(5)+yl(6)+yl(7)+yl(8))/2.;

    % Shape Functions
    ns(9) =(1. - xi.^2)*(1 - eta.^2);

    ns(5) =(1. - xi.^2)*(1. - eta   )/2. - ns(9)/2.;
    ns(6) =(1. + xi   )*(1. - eta.^2)/2. - ns(9)/2.;
    ns(7) =(1. - xi.^2)*(1. + eta   )/2. - ns(9)/2.;
    ns(8) =(1. - xi   )*(1. - eta.^2)/2. - ns(9)/2.;

    ns(1) =(1. - xi)*(1. - eta)/4. -(ns(5) + ns(8))/2. - ns(9)/4.;
    ns(2) =(1. + xi)*(1. - eta)/4. -(ns(5) + ns(6))/2. - ns(9)/4.;
    ns(3) =(1. + xi)*(1. + eta)/4. -(ns(6) + ns(7))/2. - ns(9)/4.;
    ns(4) =(1. - xi)*(1. + eta)/4. -(ns(7) + ns(8))/2. - ns(9)/4.;

    % Derivative of shape functions with respect to xi
    nxi(9) = -2.*xi*(1. - eta.^2);

    nxi(5) = -xi*(1. - eta)         - nxi(9)/2.;
    nxi(6) =(1. - eta.^2)/2. - nxi(9)/2.;
    nxi(7) = -xi*(1. + eta)         - nxi(9)/2.;
    nxi(8) =    -(1. - eta.^2)/2. - nxi(9)/2.;

    nxi(1) = -(1. - eta)/4. -(nxi(5) + nxi(8))/2. - nxi(9)/4.;
    nxi(2) =(1. - eta)/4. -(nxi(5) + nxi(6))/2. - nxi(9)/4.;
    nxi(3) =(1. + eta)/4. -(nxi(6) + nxi(7))/2. - nxi(9)/4.;
    nxi(4) = -(1. + eta)/4. -(nxi(7) + nxi(8))/2. - nxi(9)/4.;

    % Derivative of shape functions with respect to eta
    net(9) = -2.*eta*(1. - xi.^2);

    net(5) = -(1. - xi.^2)/2. - net(9)/2.;
    net(6) = -(1. + xi)*eta     - net(9)/2.;
    net(7) =(1. - xi.^2)/2. - net(9)/2.;
    net(8) = -(1. - xi)*eta     - net(9)/2.;

    net(1) = -(1. - xi)/4. -(net(5) + net(8))/2. - net(9)/4.;
    net(2) = -(1. + xi)/4. -(net(5) + net(6))/2. - net(9)/4.;
    net(3) =(1. + xi)/4. -(net(6) + net(7))/2. - net(9)/4.;
    net(4) =(1. - xi)/4. -(net(7) + net(8))/2. - net(9)/4.;

    % Find jacobian jac and its determinant.  Replace jac by its inverse
    jac = zeros(2,2);
    for i=1:9
        jac(1,1) = jac(1,1) + nxi(i)*xl(i);
        jac(1,2) = jac(1,2) + nxi(i)*yl(i);
        jac(2,1) = jac(2,1) + net(i)*xl(i);
        jac(2,2) = jac(2,2) + net(i)*yl(i);
    end;
    detjac = jac(1,1)*jac(2,2) - jac(2,1)*jac(1,2);
    tmp = jac(1,1)/detjac;
    jac(1,1) =  jac(2,2)/detjac;
    jac(1,2) = -jac(1,2)/detjac;
    jac(2,1) = -jac(2,1)/detjac;
    jac(2,2) =  tmp;

    % Form strain-displacement matrix b (zero entries are already set)
    b   = zeros(3,18);
    for j = 1: 9
        L = 2*j;
        k = L - 1;
        b(1,k) = jac(1,1)*nxi(j) + jac(1,2)*net(j);
        b(2,L) = jac(2,1)*nxi(j) + jac(2,2)*net(j);
        b(3,k) = b(2,L);
        b(3,L) = b(1,k);
    end;

    % Interpolate initial strains and thickness from corner values
    eth = zeros(1,3);
    thk = 0.;
    for i = 1: 9
        eth(1) = eth(1) + ns(i)*exi(i);
        eth(2) = eth(2) + ns(i)*eyi(i);
        eth(3) = eth(3) + ns(i)*esi(i);
        thk = thk + ns(i)*thc(i);
    end;

    return;
end % function q8_shape


function rbar_elem(n)

    fedata;
    tmp = 0.;
    for i = 1: N_UNKNOWN - NUM_RBAR*3;    
        tmp = tmp + S(i,i);
    end;
    scalef = tmp/N_UNKNOWN*0.01;
    
    n1 = ENODE(1,n);
    n2 = ENODE(2,n);
    
    x1 = X(n1);
    x2 = X(n2);
    
    y1 = Y(n1);
    y2 = Y(n2);

    y12 = y1 - y2;
    x21 = x2- x1;
    


    % Initialize the element stiffness matrix
    SE=[ 0  0   0   0   0   0   1   0   0;  ...
         0  0   0   0   0   0   0   1   0;  ... 
         0  0   0   0   0   0   y12 x21 1;  ...
         0  0   0   0   0   0   -1  0   0; ...
         0  0   0   0   0   0   0   -1  0; ...
         0  0   0   0   0   0   0   0   -1; ...
         1  0   y12 -1  0   0   0   0   0; ...
         0  1   x21 0   -1  0   0   0   0; ...
         0  0   1   0   0   -1  0   0   0].* scalef;
     
     RE(9, 1) = zeros();
     return;
end



function rod_stress(n)
    fedata;

    n1 = ENODE(1,n);
    n2 = ENODE(2,n);
    L = sqrt((X(n2)-X(n1)).^2 +(Y(n2)-Y(n1)).^2 );
    sine =(Y(n2)-Y(n1))/L;
    cosine =(X(n2)-X(n1))/L;
    strain =((U(n2)-U(n1))*cosine+(V(n2)-V(n1))*sine)/L;
    idp = P_ID(n);
    idm = M_ID(P_ID(n));
    if(TEMP_INPUT)
        delta_t = TEMPD - TREF(idm);
    else
        delta_t = 0.;
    end
    STRESS(1,n) = E_MODULUS(idm)*(strain-ALPHA(idm)*delta_t);
    ELEM_FORCE(1,n) = AREA(idp)*STRESS(1,n);
    return;
end % function rod_stress

function bar_stress(n)
    fedata;
    
    bar_elem(n);
    idp = P_ID(n);
    idm = M_ID(idp);
    
    %element ID and dimensions
    n1 = ENODE(1,n);
    n2 = ENODE(2,n);
    L = sqrt((X(n2)-X(n1)).^2 +(Y(n2)-Y(n1)).^2 );
    sine =(Y(n2)-Y(n1))/L;
    cosine =(X(n2)-X(n1))/L;
    
    %elemtnt properties
    area = AREA(idp);
    mom_iner = MOM_INERTIA(idp);
    g = G_MODULUS(idm);
    e = E_MODULUS(idm);
    h1 = H1(idp);
    h2 = H2(idp);
    
    %Load displacments
    d = [U(n1)
         V(n1)
         THETA_Z(n1)
         U(n2)
         V(n2)
         THETA_Z(n2)];
     
    
     %calculate force
     ELEM_FORCE(:,n) = SE*d-RE;
    
     f_x = -ELEM_FORCE(1,n)*cosine - ELEM_FORCE(2,n)*sine;
     
     %tensile stress
     sigma_m = f_x/area;
     
     %add bending stress at each node
     STRESS(1,n) = sigma_m + ELEM_FORCE(3,n)*h1/mom_iner;
     STRESS(2,n) = sigma_m + ELEM_FORCE(3,n)*h2/mom_iner;
     
      f_x = ELEM_FORCE(4,n)*cosine + ELEM_FORCE(5,n)*sine;
     
     %tensile stress
     sigma_m = f_x/area;
     
     STRESS(3,n) = sigma_m - ELEM_FORCE(6,n)*h1/mom_iner;
     STRESS(4,n) = sigma_m - ELEM_FORCE(6,n)*h2/mom_iner;
     
    return;
end % function bar_stress

function tria_stress(n)
    fedata;
    
    tria_elem(n);
    
        fedata;
    idp = P_ID(n);
    idm = M_ID(idp);
    if(TEMP_INPUT)
    delta_t = TEMPD - TREF(idm);
    else
    delta_t = 0.;
    end
    n1 = ENODE(1,n);
    n2 = ENODE(2,n);
    n3 = ENODE(3,n);
    x1 = X(n1);
    x2 = X(n2);
    x3 = X(n3);
    y1 = Y(n1);
    y2 = Y(n2);
    y3 = Y(n3);
    % calculate 2(area of triangle) and the area with absolute value
    a2 =(x2 - x1)*(y3 - y1) -(x3 - x1)*(y2 - y1);
    a = abs(a2)/2.;
    vol = a*THICK(idp);
    % Form [B]
    b=[(y2-y3)/a2     0.     (y3-y1)/a2     0.     (y1-y2)/a2    0.     ;  ...
    0.    (x3-x2)/a2      0.    (x1-x3)/a2     0.     (x2-x1)/a2 ;  ...
    (x3-x2)/a2 (y2-y3)/a2 (x1-x3)/a2 (y3-y1)/a2 (x2-x1)/a2 (y1-y2)/a2 ];
    % E matrix for plane stress, see Cook, page 21
    ee = E_MODULUS(idm)/(1. - PRATIO(idm)^2);
    en = ee*PRATIO(idm);
    G  = E_MODULUS(idm)/(2.*(1. + PRATIO(idm)));
    em = [ ee   en   0.  ;   ...
    en   ee   0.  ;   ...
    0.   0.   G  ] ; 
    
    
    
    
    idp = P_ID(n);
    idm = M_ID(idp);
    
    %element ID and dimensions
    n1 = ENODE(1,n);
    n2 = ENODE(2,n);
    n3 = ENODE(3,n);
   
    %elemtnt properties
    g = G_MODULUS(idm);
    e = E_MODULUS(idm);
    h1 = H1(idp);
    h2 = H2(idp);
    
    %Load displacments
    d = [U(n1)
         V(n1)
         U(n2)
         V(n2)
         U(n3)
         V(n3)];
     ep_not = [ALPHA(idm)*delta_t
                ALPHA(idm)*delta_t
                0];
      
                
    sigma = em*(b*d-ep_not);
    
    STRESS(1, n)  = sigma(1);
    STRESS(2,n) = sigma(2);
    STRESS(3,n) = sigma(3);
    
     
    return;
end % function tria_stress

function q4_stress(n)
    fedata;

    place  = [ 0.,    -0.577350269189626  -0.774596669241483; ...
               0.,     0.577350269189626,  0.; ...
               0.,     0.,                 0.774596669241483 ];
    wgt =  [ 2.,      1.,     0.555555555555556; ...
             0.,      1.,     0.888888888888889; ... 
             0.,      0.,     0.555555555555556 ];

    % define local variables
    % em(i,j) = modulus matrix
    % thc(i) = nodal thicknesses
    % b(i,j) = strain displacement matrix
    % bte(i,j) = transpose of b(i,j) times em(i,j)
    % place(i,j) = sampling points for gaussian quadrature
    %              j = order of quadrature to be used (nguass=1,2,3)
    %              i = 1, 2, or 3 for the sampling points
    % xi = a particular gaussian quadrature sampling point
    % eta = a particular gaussian quadrature sampling point
    % wgt(i,j) = wieghts for gaussian quadrature
    %              j = order of quadrature to be used (nguass=1,2,3)
    %              i = 1, 2, or 3 for the weight
    % n(i) = shape function ni, i=1, 2, 3, or 4
    %        n(i) is computed by subroutine shape at specified values
    %        of xi and eta
    % exi(i) = initial strain in the x direction at node i, i=1,2,3,4
    % eyi(i) = initial strain in the y direction at node i, i=1,2,3,4
    % esi(i) = initial shear strain at node i, i=1,2,3,4
    % eth(i) = interpolated strains at a gauss point
    % thc(i) = thickness of element at node i, i=1,2,3,4
    % thk = thickness of element at specified values of xi and eta

    % number of nodes in the element
    nn = 4;
    % number of degrees of freedom in the element
    nd = nn*2;
    idp = P_ID(n);
    idm = M_ID(idp);
    if(TEMP_INPUT)
        delta_t = TEMPD - TREF(idm);
    else
        delta_t = 0.;
    end
    
    
    %element ID and dimensions
    n1 = ENODE(1,n);
    n2 = ENODE(2,n);
    n3 = ENODE(3,n);
    n4 = ENODE(4,n);
    
    % store node numbers and their x and y coordinates
    for i = 1: 4
        xl(i)  = X(ENODE(i,n));
        yl(i)  = Y(ENODE(i,n));
    end

    % initialize the nodal thickness and initial strains
    for i=1:4
        thc(i) = THICK(idp);
        exi(i) = ALPHA(idm)*delta_t;
        eyi(i) = ALPHA(idm)*delta_t;
        esi(i) = 0.;
    end

    % initialize the modulus matrix
    % E matrix for plane stress, see Cook, page 21
    ee = E_MODULUS(idm)/(1. - PRATIO(idm)^2);
    en = ee*PRATIO(idm);
    G  = E_MODULUS(idm)/(2.*(1. + PRATIO(idm)));
    em = [ ee   en   0.  ;   ...
           en   ee   0.  ;   ...
           0.   0.   G  ] ;

    %Load Displacements
    d = [U(n1)
         V(n1)
         U(n2)
         V(n2)
         U(n3)
         V(n3)
         U(n4)
         V(n4)];
     ep_not = [ALPHA(idm)*delta_t
                ALPHA(idm)*delta_t
              0];
    %Calculate Damage at Guass Points      
    xi = 1/sqrt(3);
    eta = 1/sqrt(3);
    [eth,thk,detjac,ns,b]=q4_shape(xi,eta,xl,yl,exi,eyi,esi,thc);
    S3 = em*(b*d-ep_not);
     
    xi = 1/sqrt(3);
    eta = -1/sqrt(3);
   [eth,thk,detjac,ns,b]=q4_shape(xi,eta,xl,yl,exi,eyi,esi,thc);
    S2 = em*(b*d-ep_not);
    
     xi = -1/sqrt(3);
    eta = 1/sqrt(3);
   [eth,thk,detjac,ns,b]=q4_shape(xi,eta,xl,yl,exi,eyi,esi,thc);
    S4 = em*(b*d-ep_not);
    
     xi = -1/sqrt(3);
    eta = -1/sqrt(3);
   [eth,thk,detjac,ns,b]=q4_shape(xi,eta,xl,yl,exi,eyi,esi,thc);
    S1 = em*(b*d-ep_not);
    
     xi = 0;
    eta = 0;
   [eth,thk,detjac,ns,b]=q4_shape(xi,eta,xl,yl,exi,eyi,esi,thc);
    S_center = em*(b*d-ep_not);
    
    
    %Interpolate to nodes using shape function
    r = -sqrt(3);
    s = -sqrt(3);
    n1 =(1. - r)*(1. - s)/4.;
    n2 =(1. + r)*(1. - s)/4.;
    n3 =(1. + r)*(1. + s)/4.;
    n4 =(1. - r)*(1. + s)/4.;
    
    SN1 = n1.*S1+n2.*S2+n3.*S3+n4.*S4;
    
    r = -sqrt(3);
    s = -sqrt(3);
    n1 =(1. - r)*(1. - s)/4.;
    n2 =(1. + r)*(1. - s)/4.;
    n3 =(1. + r)*(1. + s)/4.;
    n4 =(1. - r)*(1. + s)/4.;
    
    SN1 = n1.*S1+n2.*S2+n3.*S3+n4.*S4;
    
    r = sqrt(3);
    s = -sqrt(3);
    n1 =(1. - r)*(1. - s)/4.;
    n2 =(1. + r)*(1. - s)/4.;
    n3 =(1. + r)*(1. + s)/4.;
    n4 =(1. - r)*(1. + s)/4.;
    
    SN2 = n1.*S1+n2.*S2+n3.*S3+n4.*S4;
    
    r = sqrt(3);
    s = sqrt(3);
    n1 =(1. - r)*(1. - s)/4.;
    n2 =(1. + r)*(1. - s)/4.;
    n3 =(1. + r)*(1. + s)/4.;
    n4 =(1. - r)*(1. + s)/4.;
    
    SN3 = n1.*S1+n2.*S2+n3.*S3+n4.*S4;
    
    r = -sqrt(3);
    s = sqrt(3);
    n1 =(1. - r)*(1. - s)/4.;
    n2 =(1. + r)*(1. - s)/4.;
    n3 =(1. + r)*(1. + s)/4.;
    n4 =(1. - r)*(1. + s)/4.;
    
    SN4 = n1.*S1+n2.*S2+n3.*S3+n4.*S4;
    
    
    %store stresss values
    STRESS(1,n) = SN1(1);
    STRESS(2,n) = SN1(2);
    STRESS(3,n) = SN1(3);
    STRESS(4,n) = SN2(1);
    STRESS(5,n) = SN2(2);
    STRESS(6,n) = SN2(3);
    STRESS(7, n) = SN3(1);
    STRESS(8, n) = SN3(2);
    STRESS(9, n) = SN3(3);
    STRESS(10, n) = SN4(1);
    STRESS(11, n) = SN4(2);
    STRESS(12, n) = SN4(3);
    STRESS(13, n) = S_center(1);
    STRESS(14, n) = S_center(2);
    STRESS(15, n) = S_center(3);

    return;

end

function qm6_stress(n)
    fedata;

    % We need to recover the condensed d.o.f. values.
    % generate the element stiffness matrix and load vector
    qm6_elem(n);
    % Solve [kcc]{dc} = - [kcr]{dr} + {rc}
    % Note that per Cook, 'Concepts and Applications of Finite Element
    % Analysis', 3rd Edition, page 236, calculated stresses are more
    % accurate if {rc} is set to zero during recovery.
    RE(9)  = 0.;
    RE(10) = 0.;
    RE(11) = 0.;
    RE(12) = 0.;
    % Follow the Algorithm in Cook, 3rd Ed. page 231
    de = zeros(12,1);
    de(1)  = U(ENODE(1,n));
    de(2)  = V(ENODE(1,n));
    de(3)  = U(ENODE(2,n));
    de(4)  = V(ENODE(2,n));
    de(5)  = U(ENODE(3,n));
    de(6)  = V(ENODE(3,n));
    de(7)  = U(ENODE(4,n));
    de(8)  = V(ENODE(4,n));
    de(9)  = 0.;  % not known yet
    de(10) = 0.;  % not known yet
    de(11) = 0.;  % not known yet
    de(12) = 0.;  % not known yet
    for j = 1: 4
        i = 8 + j;
        tmp = SE(i,1:i-1)*de(1:i-1);
        de(i) =(RE(i) - tmp)/SE(i,i);
    end

    % set of terms to allow evaluation of [B]
    for i = 1: 4
        xl(i)  = X(ENODE(i,n));
        yl(i)  = Y(ENODE(i,n));
    end
    idp = P_ID(n);
    idm = M_ID(idp);
    if(TEMP_INPUT)
        delta_t = TEMPD - TREF(idm);
    else
        delta_t = 0.;
    end
    % initialize the nodal thickness and initial strains
    for i=1:4
        thc(i) = THICK(idp);
        exi(i) = ALPHA(idm)*delta_t;
        eyi(i) = ALPHA(idm)*delta_t;
        esi(i) = 0.;
    end
    eth = zeros(3,1);

    % initialize the modulus matrix
    % E matrix for plane stress, see Cook, page 21
    ee = E_MODULUS(idm)/(1. - PRATIO(idm)^2);
    en = ee*PRATIO(idm);
    G  = E_MODULUS(idm)/(2.*(1. + PRATIO(idm)));
    em = [ ee   en   0.  ;   ...
           en   ee   0.  ;   ...
           0.   0.   G  ] ;

    sr3 = 1./sqrt(3.);
    for k = 1: 5
        if(k == 1)
            xi = -sr3;
            eta = -sr3;
        elseif(k == 2) 
            xi =  sr3;
            eta = -sr3;
        elseif(k == 3) 
            xi =  sr3;
            eta =  sr3;
        elseif(k == 4) 
            xi = -sr3;
            eta =  sr3;
        elseif(k == 5) 
            xi = 0.;
            eta =  0.;
        end
        % compute b matrix
        [eth,thk,detjac,ns,b,djz]=qm6_shape(xi,eta,xl,yl,exi,eyi,esi,thc);

        % compute strains = [B]{d} (note, zero multiplies are skipped)
        eps = zeros(3,1);
    %     eps = b(1:3,1:12)*de(1:12);
        eps = b*de;

        % compute stresses = [E]({e}-{ei})
        sxxg(k) = em(1,1:3)*(eps-eth');
        syyg(k) = em(2,1:3)*(eps-eth');
        sxyg(k) = em(3,1:3)*(eps-eth');
    end

    % now extrapolate to get the nodal stresses
    k = 0;
    sr3 = sqrt(3.);
    for k = 1: 4
        if(k == 1)
            xi = -sr3;
            eta = -sr3;
        elseif(k == 2) 
            xi =  sr3;
            eta = -sr3;
        elseif(k == 3) 
            xi =  sr3;
            eta =  sr3;
        elseif(k == 4) 
            xi = -sr3;
            eta =  sr3;
        end
        ns1 = (1.-xi)*(1.-eta)/4.;
        ns2 = (1.+xi)*(1.-eta)/4.;
        ns3 = (1.+xi)*(1.+eta)/4.;
        ns4 = (1.-xi)*(1.+eta)/4.;
        j =(k-1)*3;
        STRESS(1+j,n) = ns1*sxxg(1) + ns2*sxxg(2) + ns3*sxxg(3 )+ ns4*sxxg(4);
        STRESS(2+j,n) = ns1*syyg(1) + ns2*syyg(2) + ns3*syyg(3 )+ ns4*syyg(4);
        STRESS(3+j,n) = ns1*sxyg(1) + ns2*sxyg(2) + ns3*sxyg(3 )+ ns4*sxyg(4);
    end
    % now compute the stress at the element center
    STRESS(13,n) = sxxg(5);
    STRESS(14,n) = syyg(5);
    STRESS(15,n) = sxyg(5);

    return;
end % function qm6_stress


function q8_stress(n)
    fedata;

    % We need to recover the condensed d.o.f. values.
    % generate the element stiffness matrix and load vector
    q8_elem(n);
    % Solve [kcc]{dc} = - [kcr]{dr} + {rc}
    % Note that per Cook, 'Concepts and Applications of Finite Element
    % Analysis', 3rd Edition, page 236, calculated stresses are more
    % accurate if {rc} is set to zero during recovery.
    
    % Follow the Algorithm in Cook, 3rd Ed. page 231
    de = zeros(18,1);
    de(1)  = U(ENODE(1,n));
    de(2)  = V(ENODE(1,n));
    de(3)  = U(ENODE(2,n));
    de(4)  = V(ENODE(2,n));
    de(5)  = U(ENODE(3,n));
    de(6)  = V(ENODE(3,n));
    de(7)  = U(ENODE(4,n));
    de(8)  = V(ENODE(4,n));
    de(9)  = U(ENODE(5,n));
    de(10)  = V(ENODE(5,n));
    de(11)  = U(ENODE(6,n));
    de(12)  = V(ENODE(6,n));
    de(13)  = U(ENODE(7,n));
    de(14)  = V(ENODE(7,n));
    de(15)  = U(ENODE(8,n));
    de(16)  = V(ENODE(8,n));
    de(17)  = 0;
    de(18)  = 0;


    
    for j = 1:2;
    i = 16 + j;
    tmp = SE(i,1:i-1)*de(1:i-1);
    de(i) =(RE(i) - tmp)/SE(i,i);
    end

    % set of terms to allow evaluation of [B]
    for i = 1: 8
        xl(i)  = X(ENODE(i,n));
        yl(i)  = Y(ENODE(i,n));
    end
    idp = P_ID(n);
    idm = M_ID(idp);
    if(TEMP_INPUT)
        delta_t = TEMPD - TREF(idm);
    else
        delta_t = 0.;
    end
    % initialize the nodal thickness and initial strains
    for i=1:9
        thc(i) = THICK(idp);
        exi(i) = ALPHA(idm)*delta_t;
        eyi(i) = ALPHA(idm)*delta_t;
        esi(i) = 0.;
    end
    eth = zeros(3,1);

    % initialize the modulus matrix
    % E matrix for plane stress, see Cook, page 21
    ee = E_MODULUS(idm)/(1. - PRATIO(idm)^2);
    en = ee*PRATIO(idm);
    G  = E_MODULUS(idm)/(2.*(1. + PRATIO(idm)));
    em = [ ee   en   0.  ;   ...
           en   ee   0.  ;   ...
           0.   0.   G  ] ;

    sr3 = 1./sqrt(3.);
    for k = 1: 5
        if(k == 1)
            xi = -sr3;
            eta = -sr3;
        elseif(k == 2) 
            xi =  sr3;
            eta = -sr3;
        elseif(k == 3) 
            xi =  sr3;
            eta =  sr3;
        elseif(k == 4) 
            xi = -sr3;
            eta =  sr3;
        elseif(k == 5) 
            xi = 0.;
            eta =  0.;
        end
        % compute b matrix
        [eth,thk,detjac,ns,b]=q8_shape(xi,eta,xl,yl,exi,eyi,esi,thc);

        % compute strains = [B]{d} (note, zero multiplies are skipped)
        eps = zeros(3,1);
    %     eps = b(1:3,1:12)*de(1:12);
        eps = b*de;

        % compute stresses = [E]({e}-{ei})
        sxxg(k) = em(1,1:3)*(eps-eth');
        syyg(k) = em(2,1:3)*(eps-eth');
        sxyg(k) = em(3,1:3)*(eps-eth');
    end

    % now extrapolate to get the nodal stresses
    k = 0;
    sr3 = sqrt(3.);
    for k = 1: 5
        if(k == 1)
            xi = -sr3;
            eta = -sr3;
        elseif(k == 2) 
            xi =  sr3;
            eta = -sr3;
        elseif(k == 3) 
            xi =  sr3;
            eta =  sr3;
        elseif(k == 4) 
            xi = -sr3;
            eta =  sr3;
        elseif(k ==5)
            xi = 0;
            eta = 0;
        end
        

        ns1 = (1.-xi)*(1.-eta)/4.;
        ns2 = (1.+xi)*(1.-eta)/4.;
        ns3 = (1.+xi)*(1.+eta)/4.;
        ns4 = (1.-xi)*(1.+eta)/4.;
        j =(k-1)*3;
        STRESS(1+j,n) = ns1*sxxg(1) + ns2*sxxg(2) + ns3*sxxg(3 )+ ns4*sxxg(4);
        STRESS(2+j,n) = ns1*syyg(1) + ns2*syyg(2) + ns3*syyg(3 )+ ns4*syyg(4);
        STRESS(3+j,n) = ns1*sxyg(1) + ns2*sxyg(2) + ns3*sxyg(3 )+ ns4*sxyg(4);
    end

    return;
end % function q8_stress

function disp_output()
    fedata;

    fprintf(OUTFID,'%s \n', 'MAE-6010 Class; Finite Element Analysis of:');
    fprintf(OUTFID,'%s \n \n', TITLE);
    fprintf(OUTFID,'%s \n', 'Input Data Summary');
    fprintf(OUTFID,'Number of Nodes    = %8d \n', NUM_NODE);
    fprintf(OUTFID,'Number of Elements = %8d \n', NUM_ELEM);
    fprintf(OUTFID,'Constrained Degrees of Freedom \n');
    fprintf(OUTFID,'   Node Number     Direction \n');
    for i = 1: NUM_SPC
        fprintf(OUTFID,'%9d %14d \n', SPC_NODE(i),  SPC_DIR(i));
    end
    fprintf(OUTFID,'%s %5d \n', 'Number of Loaded Nodes =', NUM_LOADS);
    fprintf(OUTFID,'%s \n', '    Node ID   X Force        Y Force');
    for i = 1: NUM_LOADS
        fprintf(OUTFID,'%8d %12.5g %14.5g \n',LOADED_NODES(i),F_X(i),F_Y(i));
    end
    fprintf(OUTFID,'Accel. due to Gravity in the X = %8.5g \n',GX);
    fprintf(OUTFID,'Accel. due to Gravity in the Y = %8.5g \n',GY);
    if(TEMP_INPUT)
        fprintf(OUTFID,'%s %13.5g \n','Default Temperature =',TEMPD);
    else
        fprintf(OUTFID,'%s \n','No default temperature input');
    end

    fprintf(OUTFID,'\n %s \n','Node Displacements');
    fprintf(OUTFID,'%s \n','  Node');
    fprintf(OUTFID,'%s \n','   ID          u            v        theta Z');
    for i = 1: NUM_NODE
        fprintf(OUTFID,'%5d  %11.5g  %11.5g  %11.5g \n',i,U(i),V(i),THETA_Z(i));
    end
    disp('*****');
    return;
end % function disp_output







function rod_elem_out()
    fedata;

    found = false;
    for i = 1: NUM_ELEM
        if(ELEM_TYPE(i) == 1)
            found = true;
            break;
        end
    end
    % if none present, skip printing anything
    if(~found)
        return;
    end
    fprintf(OUTFID,'\nForces and Stresses in ROD Elements \n');
    fprintf(OUTFID,' Elem. \n');
    fprintf(OUTFID,'  ID       Force       Stress \n');
    for i = 1: NUM_ELEM
        ELEM_TYPE(i);
        if(ELEM_TYPE(i) == 1)
            fprintf(OUTFID,'  %3d  %10.5g  %10.5g \n',i,ELEM_FORCE(1,i),STRESS(1,i));
        end
    end

    return;
end % function rod_elem_out




function bar_elem_out()
    fedata;

    found = false;
    for n = 1: NUM_ELEM
        if(ELEM_TYPE(n) == 2)
            found = true;
            break;
        end
    end
    if(~found)
        return; % if none present, skip printing
    end

    fprintf(OUTFID,'\nForces at the Ends of Bar Elements in the Global Coordinate System \n');
    fprintf(OUTFID,' Elem. Node \n');
    fprintf(OUTFID,'  ID    ID        Fx          Fy          Mz \n');
    for n = 1: NUM_ELEM
        if(ELEM_TYPE(n) == 2)
            fprintf(OUTFID,'%4d  %4d %11.5g %11.5g %11.5g \n',n,ENODE(1,n),ELEM_FORCE(1:3,n));
            fprintf(OUTFID,'%4d  %4d %11.5g %11.5g %11.5g \n',n,ENODE(2,n),ELEM_FORCE(4:6,n));
        end
    end
    fprintf(OUTFID,'\nStresses at the Ends of Bar Elements in the Element Coordinate System \n');
    fprintf(OUTFID,'at Distances H1 and H2 from the Neutral Axis \n');
    fprintf(OUTFID,' Elem. Node \n');
    fprintf(OUTFID,'  ID   ID        H1           H2 \n');
    for n = 1: NUM_ELEM
        if(ELEM_TYPE(n) == 2)
            fprintf(OUTFID,'%4d %4d %11.5g  %11.5g \n',n,ENODE(1,n),STRESS(1:2,n)); 
            fprintf(OUTFID,'%4d %4d %11.5g  %11.5g \n',n,ENODE(2,n),STRESS(3:4,n));
        end
    end
    return;
end % function bar_elem_out









function tria_elem_out()
    fedata;

    found = false;
    for i = 1: NUM_ELEM
        if(ELEM_TYPE(i) == 3)
            found = true;
            break;
        end
    end
    % if none present, skip printing
    if(~found)
        return;
    end
    fprintf(OUTFID,'\nStresses in Constant Strain Triangle Elements ');
    fprintf(OUTFID,'in the Global Coordinate System \n');
    fprintf(OUTFID,' Elem.     Sigma       Sigma       Sigma       Sigma       Sigma    Von Mises\n');
    fprintf(OUTFID,'  ID        xx          yy          xy           1           2        Stress\n');
    for n = 1: NUM_ELEM
        if(ELEM_TYPE(n) == 3)
            % Calculate Max. Princpal Stresses & angle
            sx  = STRESS(1,n);
            sy  = STRESS(2,n);
            sxy = STRESS(3,n);
            tau_max = sqrt(((sx - sy)/2.).^2 + sxy.^2);
            s1 =(sx + sy)/2. + tau_max;
            s2 =(sx + sy)/2. - tau_max;
            svm = sqrt((s1-s2).^2 + s1.^2 + s2.^2)/sqrt(2.);
            fprintf(OUTFID,'%4d %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g\n',n,sx,sy,sxy,s1,s2,svm);
        end
    end
    return;
end % function tria_elem_out









function q4_elem_out()
fedata;

found = false;
for i = 1: NUM_ELEM
    if(ELEM_TYPE(i) == 4)
        found = true;
        break;
    end
end
% if none present, skip printing
if(~found)
    return;
end
if(USE_QM6)
    fprintf(OUTFID,'\nStresses in Quadrilateral QM6 Elements ');
else
    fprintf(OUTFID,'\nStresses in Quadrilateral Q4 Elements ');
end
fprintf(OUTFID,'in the Global Coordinate System \n');
fprintf(OUTFID,'Elem.           Sigma      Sigma      Sigma      Sigma      Sigma   Von Mises \n');
fprintf(OUTFID,' ID  node        xx         yy         xy          1          2      Stress \n');
for n = 1: NUM_ELEM
    if(ELEM_TYPE(n) == 4)
        k = 0;
        for i = 1: 5
            if(i < 5)
                str = num2str(ENODE(i,n));
                L = length(str);
                str = [blanks(3-L),str,blanks(2)];
            else
                str = 'Center';
            end
            % Calculate Max. Princpal Stresses & angle
            sx  = STRESS(1+k,n);
            sy  = STRESS(2+k,n);
            sxy = STRESS(3+k,n);
            tau_max = sqrt(((sx - sy)/2.).^2 + sxy.^2);
            s1 =(sx + sy)/2. + tau_max;
            s2 =(sx + sy)/2. - tau_max;
            svm = sqrt((s1-s2).^2 + s1.^2 + s2.^2)/sqrt(2.);
            fprintf(OUTFID,'%3d %6s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n', ...
                            n,str,sx,sy,sxy,s1,s2,svm);
            k = k + 3;
        end; i = 5+1;
    end
end
% Now compute averaged nodal stresses
fprintf(OUTFID,'\nAveraged Nodal Stresses in Q4 or QM6 Elements ');
fprintf(OUTFID,'in the Global Coordinate System \n');
fprintf(OUTFID,'Node     Sigma      Sigma      Sigma      Sigma      Sigma   Von Mises \n');
fprintf(OUTFID,' ID       xx         yy         xy          1          2      Stress \n');
for i = 1: NUM_NODE
    sx  = 0.;
    sy  = 0.;
    sxy = 0.;
    cnt = 0;
    for k = 1:NUM_ELEM
        if (ELEM_TYPE(k) == 4)
            for j = 1:4
               if (i == ENODE(j,k) ) 
                   cnt = cnt + 1;
                   sx = sx + STRESS(3*j-2,k);
                   sy = sy + STRESS(3*j-1,k);
                   sxy = sxy + STRESS(3*j,k);
               end
            end
        end
    end
    sx = sx/cnt;
    sy = sy/cnt;
    sxy = sxy/cnt;
    tau_max = sqrt(((sx - sy)/2.).^2 + sxy.^2);
    s1 =(sx + sy)/2. + tau_max;
    s2 =(sx + sy)/2. - tau_max;
    svm = sqrt((s1-s2).^2 + s1.^2 + s2.^2)/sqrt(2.);
    fprintf(OUTFID,'%3d %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n', ...
                    i,sx,sy,sxy,s1,s2,svm);
end
return;
end % function Q4_elem_out









function q8_elem_out()
fedata;

found = false;
for i = 1: NUM_ELEM
    if(ELEM_TYPE(i) == 5)
        found = true;
        break;
    end
end
% if none present, skip printing
if(~found)
    return;
end
fprintf(OUTFID,'\nStresses in Quadratic Quadrilateral (Q8) Elements \n');
fprintf(OUTFID,'in the Global Coordinate System \n');
fprintf(OUTFID,'Elem.           Sigma      Sigma      Sigma      Sigma      Sigma   Von Mises \n');
fprintf(OUTFID,' ID  node        xx         yy         xy          1          2      Stress \n');
for n = 1: NUM_ELEM
    if(ELEM_TYPE(n) == 5)
        k = 0;
        for i = 1: 5
            if(i < 5)
                str = num2str(ENODE(i,n));
                L = length(str);
                str = [blanks(3-L),str,blanks(2)];
            else
                str = 'Center';
            end
            % Calculate Max. Princpal Stresses & angle
            sx  = STRESS(1+k,n);
            sy  = STRESS(2+k,n);
            sxy = STRESS(3+k,n);
            tau_max = sqrt(((sx - sy)/2.).^2 + sxy.^2);
            s1 =(sx + sy)/2. + tau_max;
            s2 =(sx + sy)/2. - tau_max;
            svm = sqrt((s1-s2).^2 + s1.^2 + s2.^2)/sqrt(2.);
            fprintf(OUTFID,'%3d %6s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n',n,str,sx,sy,sxy,s1,s2,svm);
            k = k + 3;
        end; i = 5+1;
    end
end
return;
end % function Q8_elem_out






function reactf_output()
fedata;

fprintf(OUTFID,'\nReaction Forces in the Displacement Coordinate System \n');
fprintf(OUTFID,'\n  Node \n');
fprintf(OUTFID,'   ID    Direction    Force\n');
d{1} = ' X   ';
d{2} = ' Y   ';
d{3} = 'Error';
d{4} = 'Error';
d{5} = 'Error';
d{6} = 'MZ   ';
for n = 1: NUM_NODE
    for j = 1: 6
        k = ID(j,n);
        if(k > N_UNKNOWN )
            fprintf(OUTFID,'  %3d      %s   %12.5g \n',n,d{j},REACT_FORCE(k));
        end
    end
end
return;
end % function reactf_output








function trans_element(node_id, node_cnt, nterms, etype)
fedata;

% Find the c.s. identification number
found = false;
for i = 1: NUM_CS
    if(CS_ID(i) == DCS(node_id))
        % it holds which c.s system to Use
        it = i;
        found = true;
    end
end
if(~found)
    str = sprintf('Error!  Node %d6 %s \n',node_id, ...
        ' refers to a displacement coordinate system that could not be found.');
    display(str);
    str = sprintf('%s %d5 /n','The displacement c.s. = ',DCS(node_id));
    display(str);
    ERROR_F = true;
    return;
end
t_mat = eye(nterms,nterms); % make an identity matrix
% Form the transformation matrix by assigning the rotational transformation
% terms to the identity matrix 
if(etype == 2)
    nt = 3;  % bar element, rotational dof present
else
    nt = 2;
end
% location of first term to be transformed
first = node_cnt*nt - nt + 1;
for i = 1: nt
    for j = 1: nt
        t_mat(first+i-1,first+j-1) = T_CS(i,j,it);
    end
end
SE(1:nterms,1:nterms) = transpose(t_mat)*SE(1:nterms,1:nterms)*t_mat;
return;
end %function trans_element






function make_cs_transform(a,b,c,ncs)
fedata;
% move array contents to indiviual variables
x1 = a(1);
y1 = a(2);
z1 = a(3);
x2 = b(1);
y2 = b(2);
z2 = b(3);
x3 = c(1);
y3 = c(2);
z3 = c(3);

% coordinates of a vector in the x'-z' plane
v1(1) = x3 - x1;
v1(2) = y3 - y1;
v1(3) = z3 - z1;
mag = sqrt(v1(1)^2 + v1(2)^2 + v1(3)^2 );
v1 = v1./mag; % make it a Unit vector

% coordinates of a vector in the z' direction
vz(1) = x2 - x1;
vz(2) = y2 - y1;
vz(3) = z2 - z1;
mag = sqrt(vz(1)^2 + vz(2)^2 + vz(3)^2 );
vz = vz./mag; % make it a Unit vector

% compute a unit vector in the y' direction by cross product
vy = cross(vz,v1);
mag = sqrt(vy(1)^2 + vy(2)^2 + vy(3)^2 );
vy = vy./mag; % make it a Unit vector

% compute a unit vector in the x' direction by cross product
vx = cross(vy,vz);
mag = sqrt(vx(1)^2 + vx(2)^2 + vx(3)^2 );
vx = vx./mag; % make it a Unit vector

% Form the rotational transformation matrix
% Note:  This transformation matrix transforms vectors
%        FROM the new coordinate system
%        TO the reference coordinate system
T_CS(1:3,1,ncs) = vx(1:3);
T_CS(1:3,2,ncs) = vy(1:3);
T_CS(1:3,3,ncs) = vz(1:3);
return;
end






function [n]=transform_disp(n)
fedata;

% Find the c.s. identification number
found = false;
for i = 1: NUM_CS
    if(CS_ID(i) == DCS(n))
        % it holds which c.s system to use
        it = i;
        found = true;
    end
end
if(~found)
    str = sprintf('Error!  Node %d has refers to a displacement',n);
    disp([str,' coordinate system which could not be found.']);
    str = sprintf('The displacement c.s. = %d',DCS(n));
    disp([str,'  This error occured will transforming displacements.']);
    ERROR_F = true;
    return;
end
% multiply by the transformation matrix
tmat = T_CS(1:3,1:3,it);
v1 = [U(n);V(n);0.];
v2 = tmat*v1;
U(n) = v2(1);
V(n) = v2(2);
return;
end % function transform_disp






function transform_force(n)
fedata;

% Find the c.s. identification number
found = false;
it = 0;
for i = 1: NUM_CS
    if(CS_ID(i) == DCS(n))
        % it holds which c.s system to use
        it = i;
        found = true;
    end
end
if(~found)
    str = sprinf('Error!  Node %d  has refers to a displacement',n);
    disp([str,'coordinate system which could not be found.']);
    str = sprinf('The displacement c.s. = %d ',DCS(n));
    disp([str,'This error occured will transforming displacements.']);
    ERROR_F = true;
    return;
end
% Since we only support plane problems in this program, we only have to
% transform forces in the X and Y directions.  We don't need to worry
% moments since the Z axis does not change direction.
% Multiply by the force terms by the TRANSPOSE of the transformation
% matrix.
kx = ID(1,n);
ky = ID(2,n);
if(or(kx==0, ky==0))
    disp('Error! Invalid d.o.f. number in sub. transform_force');
    ERROR_F = true;
    return;
end
tmat = T_CS(1:3,1:3,it);
v1 = [R(kx);R(ky);0.];
v2 = transpose(tmat)*v1;
R(kx) = v2(1);
R(ky) = v2(2);
return;
end % function transform_force





function count_cards()
fedata;

% Count the number of elements, nodes, loads, materials, properties.
node_cnt = 0;
elem_cnt = 0;
prop_cnt = 0;
mat_cnt  = 0;
cs_cnt   = 0;
load_cnt = 0;
rbar_cnt = 0;
spc_cnt  = 0;
strng = blanks(80);
done = false;
while(~done)
    strng = my_fgetl(INFID);
    if (strng == -1)  
        done = true; 
        break;
    end
    card = strtrim(strng(1:8)); % Renoves leading and trailing blanks
    ncc = length(card);
    card = [ card, blanks(8-ncc) ]; % Relace training blanks
    if( ncc == 0 )      % skip blank lines
        continue;
    elseif (strcmp(card(1:4),'CBAR')) % ***************** CBAR ****
        elem_cnt = elem_cnt + 1;
    elseif (strcmp(card(1:6),'CORD2R') )
        cs_cnt = cs_cnt + 1;
    elseif (strcmp(card(1:6),'CQUAD4'))
        elem_cnt = elem_cnt + 1;
    elseif (strcmp(card(1:6),'CQUAD8'))
        elem_cnt = elem_cnt + 1;
    elseif (strcmp(card(1:4),'CROD'))
        elem_cnt = elem_cnt + 1;
    elseif (strcmp(card(1:6),'CTRIA3'))
        elem_cnt = elem_cnt + 1;
    elseif (strcmp(card(1:5),'FORCE'))
        load_cnt = load_cnt + 1;
    elseif (strcmp(card(1:4),'GRID'))
        node_cnt = node_cnt + 1;
    elseif (strcmp(card(1:4),'MAT1'))
        i = str2num(strng( 9:16));
        if(i>mat_cnt) mat_cnt = i; end
    elseif (strcmp(card(1:4),'PBAR'))
        i = str2num(strng( 9:16));
        if(i>prop_cnt) prop_cnt = i; end
    elseif (strcmp(card(1:4),'PROD'))
        i = str2num(strng( 9:16));
        if(i>prop_cnt) prop_cnt = i; end
    elseif (strcmp(card(1:6),'PSHELL'))
        i = str2num(strng( 9:16));
        if(i>prop_cnt) prop_cnt = i; end
    elseif (strcmp(card(1:4),'RBAR'))
        elem_cnt = elem_cnt + 1;
        rbar_cnt = rbar_cnt + 1;
    elseif (strcmp(card(1:4),'SPC '))
        spc_cnt = spc_cnt + 6;
    elseif (strcmp(card(1:4),'SPC1'))
        spc_cnt = spc_cnt + 6;
    end
end     
% Initialize the ID array as described by Cook, et. al.
ID(1,1:node_cnt) = 0;  % x displacements assumed free
ID(2,1:node_cnt) = 0;  % y displacements assumed free
ID(3,1:node_cnt) = 1;  % z displacements are restrained
ID(4,1:node_cnt) = 1;  % no rotation about x allowed
ID(5,1:node_cnt) = 1;  % no rotation about y allowed
ID(6,1:node_cnt) = 1;  % z rotations assumed constrained, 
                       % but is freed each time a beam elemnt in used
ENODE        = int32(zeros(8,node_cnt));
ELEM_TYPE    = int32(zeros(1,elem_cnt));
LOADED_NODES = int32(zeros(1,load_cnt));
DCS          = int32(zeros(1,node_cnt));
CS_ID        = int32(zeros(1,cs_cnt));
P_ID         = int32(zeros(1,elem_cnt));
M_ID         = int32(zeros(1,prop_cnt));
PHYS_PROP    = int32(zeros(1,prop_cnt));
MAT_ID       = int32(zeros(1,mat_cnt));
SPC_NODES    = int32(zeros(1,spc_cnt));
SPC_DIR      = int32(zeros(1,spc_cnt));
RBAR_ID      = int32(zeros(1,rbar_cnt));
X            = zeros(1,node_cnt);
Y            = zeros(1,node_cnt);
U            = zeros(1,node_cnt);
V            = zeros(1,node_cnt);
THETA_Z      = zeros(1,node_cnt);
DCS          = zeros(1,node_cnt);
F_X          = zeros(1,load_cnt);
F_Y          = zeros(1,load_cnt);
AREA         = zeros(1,prop_cnt);
MOM_INERTIA  = zeros(1,prop_cnt);
H1           = zeros(1,prop_cnt);
H2           = zeros(1,prop_cnt);
KSHEAR       = zeros(1,prop_cnt);
THICK        = zeros(1,prop_cnt);
E_MODULUS    = zeros(1,mat_cnt);
G_MODULUS    = zeros(1,mat_cnt);
RHO          = zeros(1,mat_cnt);
ALPHA        = zeros(1,mat_cnt);
TREF         = zeros(1,mat_cnt);
PRATIO       = zeros(1,mat_cnt);
T_CS         = zeros(3,3,cs_cnt);
CONSTR_DISP  = zeros(1,spc_cnt);
STRESS       = zeros(16,elem_cnt);
ELEM_FORCE   = zeros(6,elem_cnt);
REACT_FORCE  = zeros(1,node_cnt);
end





function draw_mesh()
fedata;
figure;
clf;
axes('DataAspectRatioMode','manual');
hold on;
for n=1:NUM_ELEM
    typ = ELEM_TYPE(n);
    [xp,yp,npt]=get_xy(0.,n,typ);
    if(typ==1)
        plot(xp,yp,'-b.','LineWidth',2,'MarkerEdgeColor','k');
    elseif(typ==2)
        plot(xp,yp,'-m.','LineWidth',2,'MarkerEdgeColor','k');
    elseif(typ==6)
        plot(xp,yp,'-y.','LineWidth',3,'MarkerEdgeColor','k');
    else
        plot(xp,yp,'-k.','LineWidth',1,'MarkerEdgeColor','k');
        fill(xp,yp,[.6 1 .6]);
    end
    xc=0.;
    yc=0.;
    for i=1:npt
        xc=xc+xp(i);
        yc=yc+yp(i);
    end
    xc=xc/double(npt);
    yc=yc/double(npt);
    str=num2str(n);
    text(xc,yc,str,'Color','r');
    xp=zeros(1,2);
    yp=zeros(1,2);
end
for n=1:NUM_NODE
    xpp=X(n);
    ypp=Y(n);
    str=[' ',num2str(n)];
    text(xpp,ypp,str);
end
xp=zeros(1,2);
yp=zeros(1,2);
hold off
end




function draw_deformed()
fedata;
xmax=max(X);
xmin=min(X);
ymax=max(Y);
ymin=min(Y);
dc=max(xmax-xmin,ymax-ymin);
umax=max(U);
umin=min(U);
vmax=max(V);
vmin=min(V);
dd=max(umax-umin,vmax-vmin);
sf=dc/dd*.15;  % Scale factor applied to displacements
figure;
clf;
axes('DataAspectRatioMode','manual');
hold on;
for n=1:NUM_ELEM
    typ = ELEM_TYPE(n);
    [xp,yp,npt]=get_xy(sf,n,typ);
    if(typ==1)
        plot(xp,yp,'-b','LineWidth',2);
    elseif(typ==2)
        plot(xp,yp,'-m','LineWidth',2);
    elseif(typ==6)
        plot(xp,yp,'-y','LineWidth',3);
    else
        plot(xp,yp,'-k','LineWidth',1);
        fill(xp,yp,[.6 1 .6]);
    end
    xp=zeros(1,2);
    yp=zeros(1,2);
end
colormap;
hold off
end




function [xp,yp,npt]=get_xy(sf,n,typ)
fedata;
npt=typ;
if(typ == 1) npt=2; end
if(typ == 5) npt=8; end
if(typ == 6) npt=2; end
if(npt<8)
    for i=1:npt
        xp(i)=X(ENODE(i,n))+U(ENODE(i,n))*sf;
        yp(i)=Y(ENODE(i,n))+V(ENODE(i,n))*sf;
    end
else
    xp(1)=X(ENODE(1,n))+U(ENODE(1,n))*sf;
    yp(1)=Y(ENODE(1,n))+V(ENODE(1,n))*sf;
    xp(2)=X(ENODE(5,n))+U(ENODE(5,n))*sf;
    yp(2)=Y(ENODE(5,n))+V(ENODE(5,n))*sf;
    xp(3)=X(ENODE(2,n))+U(ENODE(2,n))*sf;
    yp(3)=Y(ENODE(2,n))+V(ENODE(2,n))*sf;
    xp(4)=X(ENODE(6,n))+U(ENODE(6,n))*sf;
    yp(4)=Y(ENODE(6,n))+V(ENODE(6,n))*sf;
    xp(5)=X(ENODE(3,n))+U(ENODE(3,n))*sf;
    yp(5)=Y(ENODE(3,n))+V(ENODE(3,n))*sf;
    xp(6)=X(ENODE(7,n))+U(ENODE(7,n))*sf;
    yp(6)=Y(ENODE(7,n))+V(ENODE(7,n))*sf;
    xp(7)=X(ENODE(4,n))+U(ENODE(4,n))*sf;
    yp(7)=Y(ENODE(4,n))+V(ENODE(4,n))*sf;
    xp(8)=X(ENODE(8,n))+U(ENODE(8,n))*sf;
    yp(8)=Y(ENODE(8,n))+V(ENODE(8,n))*sf;
end
end