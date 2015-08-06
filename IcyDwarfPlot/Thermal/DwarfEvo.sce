    xdel(winsid());
    
    TITLE = string("Ceres interior temperatures");
    infile = 'Thermal.txt';
    Radius_max = 500;
    NR = 200;         //Number of shells
    NT = 501;         //Number of timesteps (incl. zeroth)
    timestep = 0.01;  //Timestep of the sim, specified in IcyDwarfInput.txt
    G = 6.67e-11;
    Tmin = 3000;
    Tmax = 0;
    
    //Initialisations
    time = zeros(NT);
    Radius = zeros(NR);
    Temp = zeros(NT,NR);
    Mrock = zeros(NT,NR);
    Mh2os = zeros(NT,NR);
    Madhs = zeros(NT,NR);
    Mh2ol = zeros(NT,NR);
    Mnh3l = zeros(NT,NR);
    Kappa = zeros(NT,NR);
    Xhydr = zeros(NT,NR);

    //Read the input file (i.e., output file of the thermal code)
    fid = mopen(infile, "r");
    if (fid == -1) then
        error("Cannot open input file");
    end
    
    for (t = 1:1:NT)
        if (t>1) 
            time(t) = time(t-1) + timestep;
        end
        for (r = 1:1:NR)
             [num_read, Radius(r), Temp(t,r), Mrock(t,r), Mh2os(t,r), Madhs(t,r), Mh2ol(t,r), Mnh3l(t,r), Kappa(t,r), Xhydr(t,r)] = mfscanf(fid, "%f %f %f %f %f %f %f %f %f");
             if (Temp(t,r) < Tmin) then
                 Tmin = Temp(t,r);
             end
             if (Temp(t,r) > Tmax) then
                 Tmax = Temp(t,r);
             end
        end
    end    
    mclose(fid);
    
    MR = [0   ,0   ;
    
          0.07,0.2 ;
          0.12,0.2 ;
          
          0.139,0.9;
          0.14,0.7 ;
          
          0.19,1   ;
          0.2 ,0.2 ;
          
          0.32,0.35;
          0.38,0.8 ;
          
          0.55,0.8 ;
          0.65,0.5 ;
          
          1   ,0.5 ]; // Colormap red matrix
          
    MG = [0   ,0   ;
    
          0.07,0.2 ;
          0.12,0.75;
          
          0.139,0.9;
          0.14,0.3 ;
          
          0.19,1   ;
          0.2 ,0.3 ;
          
          0.32,0.3 ;
          0.38,0.8 ;
          
          0.55,0.8 ;
          0.65,0   ;
          
          1   ,0   ]; // Colormap green matrix
          
    MB = [0   ,0   ;
    
          0.07,0.5 ;
          0.12,0.5 ;
          
          0.139,0.9;
          0.14,0   ;
          
          0.19,1   ;
          0.2 ,0   ;
          
          0.32,0   ;
          0.38,0.3 ;
          
          0.55,0.35;
          0.65,0   ;
          
          1   ,0   ]; // Colormap blue matrix

    //Create figure
    
    f = scf();
    f.figure_name = TITLE;
    f.color_map = RGBmatrices(2048,MR,MG,MB);
    f.background = color(255,255,255);
    f.figure_size = [800,450];
    nz = [100 176 200 273 300 500 700 850 1000 1500 2000 2500];
    
    //title(TITLE,'fontsize',5);
    xlabel("Time after CAI formation (Gyr)",'fontsize',3);
    ylabel("Radius (km)",'fontsize',3);
    YL=[0:200:2000]; 
    TL=["0","200",'400","600","800","1000","1200","1400","1600","1800","2000"]; 
    colorbar(0, 2300, [1 2048], "%.0f");
    cb = gcf(); 
    cb.children(1).auto_ticks(2)="off"; 
    cb.children(1).y_ticks = tlist(["ticks","locations","labels"], YL, TL); 
    Sgrayplot(time,Radius,Temp,rect=[0,0,time(NT),Radius(NR)],zminmax=[0 2300]);
    
    
    xset('font size',0.4);
    xset("fpf","%.0f");
    //xset('thickness',1.5);
    contour2d(time,Radius,Temp,nz,rect=[0,0,time(NT),Radius(NR)]);
    
    a=get("current_axes")//get the handle of the newly created axes
    a.axes_visible="on"; // makes the axes visible
    a.data_bounds(:,1) = [0.01;10] ; // set positive bounds for X axis
    a.log_flags = "ln" ; // set X axes to logarithmic scale
    a.sub_ticks = [10 5];
    
    //Export to PNG file
    xs2png(f,"TempK.png");
    
    // Custom color map
    function mymap = RGBmatrices(N, rm, gm, bm)
        x = linspace(0,1, N);
        rv = interp1( rm(:,1), rm(:,2), x);
        gv = interp1( gm(:,1), gm(:,2), x);
        mv = interp1( bm(:,1), bm(:,2), x);
        mymap = [ rv', gv', mv'];
        // exclude invalid values that could appear
        mymap( isnan(mymap) ) = 0;
        mymap( (mymap>1) ) = 1;
        mymap( (mymap<0) ) = 0;
    endfunction
