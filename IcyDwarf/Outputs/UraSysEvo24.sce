    clear;          // Clear variables
    close;
    close;
    close;
    close;

    nmoons = 5;       // Number of moons
    NT = 451;    // Number of timesteps (incl. zeroth)
    timestep = 0.01;  // Timestep of the sim, specified in IcyDwarfInput.txt (Gyr)
    time = zeros(NT);
    time(1) = 0.01;  // Time zero of the sim (Gyr)
    tspawn = zeros(nmoons); // Time of formation (Myr)
    tspawn(1) = 10-10+1; // Miranda
    tspawn(2) = 10-10+1; // Ariel
    tspawn(3) = 10-10+1; // Umbriel
    tspawn(4) = 10-10+1;  // Titania
    tspawn(5) = 10-10+1;  // Oberon
    
    function set_line_style(thickness)
        e = gce();
        e.children.thickness = thickness;
    endfunction
    
    for (i = 1:1:nmoons)
        if (tspawn(i) > NT) 
            tspawn(i) = NT-5;
        end
    end
   
    for (m=1:1:nmoons)
        for (t = 1:1:NT) 
            if (t>1) 
                time(t) = time(t-1) + timestep;
            end
        end
    end
    
    // -----------------------------------
    // INTERNAL TEMPERATURES
    // -----------------------------------

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
    
    TITLE = string("");
    infile(1) = '0Thermal.txt';
    infile(2) = '1Thermal.txt';
    infile(3) = '2Thermal.txt';
    infile(4) = '3Thermal.txt';
    infile(5) = '4Thermal.txt';

    NR = 200;          //Number of shells
    Tmin = 5000;
    Tmax = 0;
    
    Radius = zeros(nmoons,NT,NR);
    Temp = zeros(nmoons,NT,NR);
    Mrock = zeros(nmoons,NT,NR);
    Mh2os = zeros(nmoons,NT,NR);
    Madhs = zeros(nmoons,NT,NR);
    Mh2ol = zeros(nmoons,NT,NR);
    Mnh3l = zeros(nmoons,NT,NR);
    Nu = zeros(nmoons,NT,NR);
    Famor = zeros(nmoons,NT,NR);
    Kappa = zeros(nmoons,NT,NR);
    Xhydr = zeros(nmoons,NT,NR);
    Pore = zeros(nmoons,NT,NR);
    Crack = zeros(NT,NR);
    TideHeatRate = zeros(NT,NR);
    Radiustemp = zeros(NR);
    Temptemp = zeros(NT,NR);

    //Read the thermal output files
    
    fid = zeros(nmoons);
    
    for (m = 1:1:nmoons)
        fid(m) = mopen(infile(m), "r");
        if (fid(m) == -1) then
            error("Cannot open input file");
        end
        
        for (t = 1:1:NT)
            if (t>1) 
                time(t) = time(t-1) + timestep;
            end
            for (r = 1:1:NR)
                 [num_read, Radius(m,t,r), Temp(m,t,r), Mrock(m,t,r), Mh2os(m,t,r), Madhs(m,t,r), Mh2ol(m,t,r), Mnh3l(m,t,r), Nu(m,t,r), Famor(m,t,r), Kappa(m,t,r), Xhydr(m,t,r), Pore(m,t,r), Crack(m,t,r), TideHeatRate(m,t,r)] = mfscanf(fid(m), "%f %f %f %f %f %f %f %f %f %f %f %f %f %f");
                 if (Temp(m,t,r) < Tmin) then
                     Tmin = Temp(m,t,r);
                 end
                 if (Temp(m,t,r) > Tmax) then
                     Tmax = Temp(m,t,r);
                 end
            end
        end
        
        mclose(fid(m));
    end
    
    // Set temperature color map
    MR = [0   ,0   ;
          0.02,0   ;
          0.087,0.9   ;
          0.09,0.2 ;
          0.12,0.2 ;
          
          0.1364,0.8;
          0.1366,0.7 ;
          
          0.15,1   ;
          0.2 ,0.2 ;
          
          0.32,0.35;
          0.38,0.8 ;
          
          0.55,0.8 ;
          0.65,0.5 ;
          
          1   ,0.5 ]; // Colormap red matrix
          
    MG = [0   ,0   ;
          0.02,0   ;
          0.087,0.9   ;
          0.09,0.2 ;
          0.12,0.75;
          
          0.1364,0.9;
          0.1366,0.3 ;
          
          0.15,1   ;
          0.2 ,0.3 ;
          
          0.32,0.3 ;
          0.38,0.8 ;
          
          0.55,0.8 ;
          0.65,0   ;
          
          1   ,0   ]; // Colormap green matrix
          
    MB = [0   ,0   ;
          0.02,0   ;
          0.087,0.9   ;
          0.09,0.5 ;  // Ammonia eutectic
          0.12,0.5 ;
          
          0.1364,0.8;  // Liquid water
          0.1366,0   ;
          
          0.15,1   ;
          0.2 ,0   ;
          
          0.32,0   ;  // Dehydration
          0.38,0.3 ;
          
          0.55,0.35;  // Melting
          0.65,0   ;
          
          1   ,0   ]; // Colormap blue matrix

    //Create figure
    
    f = scf();
    f.figure_name = "Temperatures";
    f.color_map = RGBmatrices(2048,MR,MG,MB);
    f.background = color(255,255,255);
    f.figure_size = [1440,400];
//    nz = [50 100 150 200 250 300 400 500 600 700 800 900 1000 1100 1200 1300];
    nz = [];
    fontsize = 3;

    for (m = 1:1:nmoons)
        for (t = 1:1:NT)
            for (r = 1:1:NR)
                Temptemp(t,r) = Temp(m,t,r);
            end
        end
        for (r = 1:1:NR)
            Radiustemp(r) = Radius(m,NT,r);
        end
        
        subplot(1,nmoons+1,m);
        if (m==1) 
            title("Miranda",'color',color("deepskyblue"),'fontsize',fontsize);
            ylabel("Radius (km)",'fontsize',fontsize);
        end
        if (m==2)
             title("Ariel",'color',color("mediumseagreen"),'fontsize',fontsize);  
        end
        if (m==3)
             title("Umbriel",'color',color("dark orange"),'fontsize',fontsize);   
        end
        if (m==4)
             title("Titania",'color',color("darkviolet"),'fontsize',fontsize);   
        end
        if (m==5)       
             title("Oberon",'color',color("deeppink"),'fontsize',fontsize); 
        end
        xlabel("Time (Gyr)",'fontsize',fontsize);
        
//        Sgrayplot(time(tspawn(m):NT),Radiustemp,Temptemp(tspawn(m):NT,:),rect=[0,0,time(50),Radiustemp(NR)],zminmax=[0 2000]);
        Sgrayplot(time(1:NT),Radiustemp,Temptemp(1:NT,:),rect=[0,0,time(10),Radiustemp(NR)],zminmax=[0 2000]);

        contour2d(time(1:NT),Radiustemp,Temptemp(1:NT,:),nz,rect=[0,0,time(10),Radiustemp(NR)]);
                                                                             
        a=get("current_axes") //get the handle of the newly created axes
        a.font_size = fontsize;
        a.axes_visible="on"; // makes the axes visible
    //    a.log_flags = "ln" ; // set X axes to logarithmic scale
        a.data_bounds(:,1) = [0;5] ; // set positive bounds for X axis
        a.sub_ticks = [8 4];
    end
    
    // Color bar, also the extra subplot helps adjust the one before
    subplot(1,nmoons+1,m+1);
    YL=[0:100:2000]; 
    TL=["0","100","200","300","400","500","600","700","800","900","1000","1100","1200","1300","1400","1500","1600","1700","1800","1900","2000"]; 
    colorbar(0, 1000, [1 1024], "%.0f");  
    cb = gcf(); 
    cb.children(1).auto_ticks(2)="off"; 
    cb.children(1).y_ticks = tlist(["ticks","locations","labels"], YL, TL); 
    
    //Export to PNG file
    xs2png(f,"TempK.png");
    
    // -----------------------------------
    // STRUCTURES
    // -----------------------------------
    
    // +2 lines to close the polygons
    Liqcore = zeros(nmoons,NT+2);
    Core = zeros(nmoons,NT+2);
    Ocean = zeros(nmoons,NT+2);
    Ice = zeros(nmoons,NT+2);
    Crust = zeros(nmoons,NT+2);
    timepoly = zeros(nmoons,NT+2);
    
    for (m = 1:1:nmoons)
        for (t = 1:1:NT)
            timepoly(m,t) = time(t);
            // Initialize as undifferentiated
            Liqcore(m,t) = 0;
            Core(m,t) = 0;
            Ocean(m,t) = 0;
            Ice(m,t) = 0;
            Crust(m,t) = Radius(m,t,NR);
            
            // Figure out outer radius of each layer
            for (r = 2:1:NR)
                if (Mrock(m,t,r) <= Mrock(m,t,r-1) & Mrock(m,t,r-1) > 0)
                    Core(m,t) = Radius(m,t,r)
                end
                if (Mrock(m,t,r) > 0 & Mrock(m,t,r-1) <= 0 & Mh2ol(m,t,r-1) <= 0)
                    Ice(m,t) = Radius(m,t,r)
                end
                if (Mh2ol(m,t,r) > 0)
                    Ocean(m,t) = Radius(m,t,r)
                end
//                if (Mrock(m,t,r) > 0 & Mh2ol(m,t,r) > 0 & Mh2ol(m,t,r-1) > 0)
//                    Liqcore(m,t) = Radius(m,t,r)
//                end
            end
            if (Mrock(m,t,NR) <= 0) 
                Ice(m,t) = Radius(m,t,NR);
            end
        end
    end
    
    // Close polygons
    for (m = 1:1:nmoons)
        timepoly(m,NT+1) = timepoly(m,NT); timepoly(m,NT+2) = min(timepoly(m,tspawn(m):NT));
        Liqcore(m,NT+1) = 0; Liqcore(m,NT+2) = 0;
        Core(m,NT+1) = 0; Core(m,NT+2) = 0;
        Ocean(m,NT+1) = 0; Ocean(m,NT+2) = 0;
        Ice(m,NT+1) = 0; Ice(m,NT+2) = 0;
        Crust(m,NT+1) = 0; Crust(m,NT+2) = 0;
    end

    //Create figure
    
    fstr = scf();
    fstr.figure_name = "Structures";
    fstr.background = color(255,255,255);
    fstr.figure_size = [1440,400];

    for (m = 1:1:nmoons)
        subplot(1,nmoons+1,m);
        if (m==1) 
            title("Miranda",'color',color("deepskyblue"),'fontsize',fontsize);
            ylabel("Radius (km)",'fontsize',fontsize);
        end
        if (m==2)
             title("Ariel",'color',color("mediumseagreen"),'fontsize',fontsize);  
        end
        if (m==3)
             title("Umbriel",'color',color("dark orange"),'fontsize',fontsize);   
        end
        if (m==4)
             title("Titania",'color',color("darkviolet"),'fontsize',fontsize);   
        end
        if (m==5)       
             title("Oberon",'color',color("deeppink"),'fontsize',fontsize); 
        end
        xlabel("Time (Gyr)",'fontsize',fontsize);
        
//        xfpoly(timepoly(m,tspawn(m):$),Crust(m,tspawn(m):$)); // ' = transpose
        xfpoly(timepoly(m,10:$),Crust(m,10:$)); // ' = transpose
        h = gce();
        h.line_mode = "off";
        h.background = color("gray50"); // Adjust color
        
        if (max(Ice(m,tspawn(m):$)) > 0) then
//            xfpoly(timepoly(m,tspawn(m):$),Ice(m,tspawn(m):$)); // ' = transpose
            xfpoly(timepoly(m,10:$),Ice(m,10:$)); // ' = transpose
            h = gce();
            h.line_mode = "off";
            h.background = color("lightblue2"); // Adjust color
        end
        
        if (max(Ocean(m,tspawn(m):$)) > 0) then
//            xfpoly(timepoly(m,tspawn(m):$),Ocean(m,tspawn(m):$)); // ' = transpose
            xfpoly(timepoly(m,10:$),Ocean(m,10:$)); // ' = transpose
            h = gce();
            h.line_mode = "off";
            h.background = color("royalblue1"); // Adjust color
        end
        
        if (max(Core(m,tspawn(m):$)) > 0) then
//            xfpoly(timepoly(m,tspawn(m):$),Core(m,tspawn(m):$)); // ' = transpose
            xfpoly(timepoly(m,10:$),Core(m,10:$)); // ' = transpose
            h = gce();
            h.line_mode = "off";
            h.background = color("chocolate"); // Adjust color
        end
        
        if (max(Liqcore(m,tspawn(m):$)) > 0) then
//            xfpoly(timepoly(m,tspawn(m):$),Liqcore(m,tspawn(m):$)); // ' = transpose
            xfpoly(timepoly(m,10:$),Liqcore(m,10:$)); // ' = transpose
            h = gce();
            h.line_mode = "off";
            h.background = color("light sea green"); // Adjust color
        end
        
        a=get("current_axes") //get the handle of the newly created axes
        a.font_size = fontsize;
        a.axes_visible="on"; // makes the axes visible
        a.log_flags = "nn" ; // lin-lin
        a.data_bounds(:,1) = [0;5] ; // set positive bounds for X axis
        a.sub_ticks = [8 4];
        
        plot2d(time,Radius(m,1:NT,NR));
        a=gca();
        poly1= a.children(1).children(1); //store polyline handle into poly1 
        poly1.foreground = color("white"); // another way to change the style
    end
    
    // Placeholder in order to scale the last subplot
    subplot(1,nmoons+1,m+1);
    plot2d(time,Radius(m,1:NT,NR));
    a=gca();
    poly1= a.children(1).children(1); //store polyline handle into poly1 
    poly1.foreground = color("white"); // another way to change the style
    a.axes_visible="off"; // makes the axes invisible
    
    //Export to PNG file
    xs2png(fstr,"Struct.png");
    
    // -----------------------------------
    // -----------------------------------
      
    // -----------------------------------
    // ORBITS
    // -----------------------------------
    
    NT = NT+10-1
    NT = NT*10;
//    NT = 4471;
    NT = NT-100+1;
    for (m = 1:1:nmoons)
        tspawn(m) = tspawn(m)+10-1;
        tspawn(m) = tspawn(m)*10;
        tspawn(m) = tspawn(m)-100+1;
        if (tspawn(i) > NT) 
            tspawn(i) = NT-5;
        end
    end
    timestep = timestep/10;
    time = zeros(NT);
    time(1) = 0.1;
    for (m=1:1:nmoons)
        for (t = 1:1:NT) 
            if (t>1) 
                time(t) = time(t-1) + timestep;
            end
        end
    end
    
    TITLEOrb = string("");
    infile(1) = '0Orbit.txt';
    infile(2) = '1Orbit.txt';
    infile(3) = '2Orbit.txt';
    infile(4) = '3Orbit.txt';
    infile(5) = '4Orbit.txt';
    
    timeorb = zeros(NT);
    aorb = zeros(nmoons,NT);
    aorbtemp = zeros(NT);
    aosc = zeros(nmoons,NT);
    eorb = zeros(nmoons,NT);
    eorbtemp = zeros(NT);
    horb = zeros(nmoons,NT);
    korb = zeros(nmoons,NT);
    resangle = zeros(nmoons,NT);
    Wtide = zeros(nmoons,NT);
    Wtidetemp = zeros(NT);
    k2Q = zeros(nmoons,NT);
    k2Qtemp = zeros(NT);
    
    moonsize = zeros(nmoons);
    
    moonsize(1) = 236;
    moonsize(2) = 579;
    moonsize(3) = 585;
    moonsize(4) = 788;
    moonsize(5) = 761;

    //Read the orbital output files
    
    fidorb = zeros(nmoons);
    
    for (m = 1:1:nmoons)
        fidorb(m) = mopen(infile(m), "r");
        if (fidorb(m) == -1) then
            error("Cannot open input file");
        end
        
        for (t = 1:1:NT)
            [num_read, timeorb(t), aorb(m,t), aosc(m,t), eorb(m,t), horb(m,t), korb(m,t), resangle(m,t), Wtide(m,t), k2Q(m,t)] = mfscanf(fidorb(m), "%f %f %f %f %f %f %f %f %f");
        end   
        
        mclose(fidorb(m)); 
    end

    //Create figure
    
    forb = scf();
    forb.figure_name = "Orbits";
    forb.background = color(255,255,255);
    forb.figure_size = [1200,600];

    // Semi-major axes
    subplot(1,4,1);
    title("Semi-major axes",'fontsize',fontsize);
    xlabel("Time (Gyr)",'fontsize',fontsize);
    ylabel("Semi-major axis (km)",'fontsize',fontsize);
    for (m = 1:1:nmoons)
        for (t = 1:1:NT)
            aorbtemp(t) = aorb(m,t);
        end
        if (m == 1) then
            plot2d(time(tspawn(m):NT),aorbtemp(tspawn(m):NT),color("deepskyblue"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 2) then
            plot2d(time(tspawn(m):NT),aorbtemp(tspawn(m):NT),color("mediumseagreen"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 3) then
            plot2d(time(tspawn(m):NT),aorbtemp(tspawn(m):NT),color("dark orange"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 4) then
            plot2d(time(tspawn(m):NT),aorbtemp(tspawn(m):NT),color("darkviolet"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 5) then
            plot2d(time(tspawn(m):NT),aorbtemp(tspawn(m):NT),color("deeppink"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
    end
    a=get("current_axes") //get the handle of the newly created axes
    a.font_size = fontsize;
    a.axes_visible="on"; // makes the axes visible
    a.log_flags = "nn" ; // lin-lin
    a.data_bounds(:,1) = [0;5] ; // set positive bounds for X axis
    a.sub_ticks = [8 4];
    
    // Overlay current semi-major axes as filled circles
    currentaorb = zeros(nmoons);
    currentaorb(1) = 129858;
    currentaorb(2) = 190390;
    currentaorb(3) = 265982;
    currentaorb(4) = 436282;
    currentaorb(5) = 583449;
                     
    currenttime = zeros(nmoons);
    for (m = 1:1:nmoons)
        currenttime(m) = 4.56;
    end
    
    if (nmoons==5)
        scatter(currenttime,currentaorb,moonsize/10,[color("black")],"fill",".");
    end
    if (nmoons <> 5) then
        mprintf("Not plotting the current orbital positions because nmoons != 5\n");
    end

    // Eccentricities
    subplot(1,4,2); 
    title("Eccentricities",'fontsize',fontsize);
    xlabel("Time (Gyr)",'fontsize',fontsize);
    ylabel("Eccentricity",'fontsize',fontsize);
    for (m = 1:1:nmoons)
        for (t = 1:1:NT)
            eorbtemp(t) = eorb(m,t);
        end
        if (m == 1) then
            plot2d(time(tspawn(m):NT),eorbtemp(tspawn(m):NT),color("deepskyblue"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 2) then
            plot2d(time(tspawn(m):NT),eorbtemp(tspawn(m):NT),color("mediumseagreen"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 3) then
            plot2d(time(tspawn(m):NT),eorbtemp(tspawn(m):NT),color("dark orange"));//,rect=[0,1e-7,time(50),1]);
        end
        if (m == 4) then
            plot2d(time(tspawn(m):NT),eorbtemp(tspawn(m):NT),color("darkviolet"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 5) then
            plot2d(time(tspawn(m):NT),eorbtemp(tspawn(m):NT),color("deeppink"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
    end
    a=get("current_axes") //get the handle of the newly created axes
    a.font_size = fontsize;
    a.axes_visible="on"; // makes the axes visible
    a.data_bounds = [0,1e-7;5,1];
    a.log_flags = "nl" ; // lin-log
    a.sub_ticks = [8 4];
    
    // Overlay current eccentricities as filled circles
    currenteorb = zeros(nmoons);
//    currenteorb(1) = 0.0196;
//    currenteorb(2) = 0.0047;
//    currenteorb(3) = 0.0001;
//    currenteorb(4) = 0.0022;
//    currenteorb(5) = 0.0010;
    
//    if (nmoons==5)
//        scatter(currenttime,currenteorb,moonsize/10,[0,0],"fill",".");
//    end
//    if (nmoons <> 5) then
//        mprintf("Not plotting the current orbital positions because nmoons != 5\n");
//    end
    
    // Tidal dissipation
    subplot(1,4,3); 
    title("Total tidal dissipation",'fontsize',fontsize);
    xlabel("Time (Gyr)",'fontsize',fontsize);
    ylabel("Tidal dissipation (W)",'fontsize',fontsize);
    for (m = 1:1:nmoons)
        for (t = 1:1:NT)
            Wtidetemp(t) = Wtide(m,t);
        end
        if (m == 1) then
            plot2d(time(tspawn(m):NT),Wtidetemp(tspawn(m):NT),color("deepskyblue"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 2) then
            plot2d(time(tspawn(m):NT),Wtidetemp(tspawn(m):NT),color("mediumseagreen"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 3) then
            plot2d(time(tspawn(m):NT),Wtidetemp(tspawn(m):NT),color("dark orange"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 4) then
            plot2d(time(tspawn(m):NT),Wtidetemp(tspawn(m):NT),color("darkviolet"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 5) then
            plot2d(time(tspawn(m):NT),Wtidetemp(tspawn(m):NT),color("deeppink"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
    end
    a=get("current_axes") //get the handle of the newly created axes
    a.font_size = fontsize;
    a.axes_visible="on"; // makes the axes visible
    a.data_bounds = [0,1e-1;5,1e12];
    a.log_flags = "nl" ; // lin-log
    a.sub_ticks = [8 4];
    
    // Equivalent k2/Q
    subplot(1,4,4); 
    title("k2/Q",'fontsize',fontsize);
    xlabel("Time (Gyr)",'fontsize',fontsize);
    ylabel("k2/Q",'fontsize',fontsize);
    for (m = 1:1:nmoons)
        for (t = 1:1:NT)
            k2Qtemp(t) = k2Q(m,t);
        end
        if (m == 1) then
            plot2d(time(tspawn(m):NT),k2Qtemp(tspawn(m):NT),color("deepskyblue"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 2) then
            plot2d(time(tspawn(m):NT),k2Qtemp(tspawn(m):NT),color("mediumseagreen"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 3) then
            plot2d(time(tspawn(m):NT),k2Qtemp(tspawn(m):NT),color("dark orange"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 4) then
            plot2d(time(tspawn(m):NT),k2Qtemp(tspawn(m):NT),color("darkviolet"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
        if (m == 5) then
            plot2d(time(tspawn(m):NT),k2Qtemp(tspawn(m):NT),color("deeppink"));//,rect=[0,1e-7,time(50),1]);
            set_line_style(2.8);
        end
    end
    a=get("current_axes") //get the handle of the newly created axes
    a.font_size = fontsize;
    a.axes_visible="on"; // makes the axes visible
    a.data_bounds = [0,1e-12;5,1e-2];
    a.log_flags = "nl" ; // lin-log
    a.sub_ticks = [8 4];
    
    //Export to PNG file
    xs2png(forb,"Orb.png");
    
    // -----------------------------------
    // RESONANCES
    // -----------------------------------
    
    MMratio = zeros(nmoons,nmoons,NT);
    
    for (m = 1:1:nmoons)
        for (n = 1:1:nmoons)
            for (t = 1:1:NT)
                MMratio(m,n,t) = (aorb(m,t)/aorb(n,t))^(-1.5);
            end
        end
    end
    
    // Overlay resonances
    res21 = zeros(5001);
    res32 = zeros(5001);
    res43 = zeros(5001);
    res54 = zeros(5001);
    res65 = zeros(5001);
    
    res31 = zeros(5001);
    res53 = zeros(5001);
    res75 = zeros(5001);
    
    res41 = zeros(5001);
    res52 = zeros(5001);
    res74 = zeros(5001);
    
    res12 = zeros(5001);
    res23 = zeros(5001);
    res34 = zeros(5001);
    res45 = zeros(5001);
    res56 = zeros(5001);
    
    res13 = zeros(5001);
    res35 = zeros(5001);
    res57 = zeros(5001);
    
    res14 = zeros(5001);
    res25 = zeros(5001);
    res47 = zeros(5001);
    
    for (i = 1:1:5001)
        res21(i) = 2/1;
        res32(i) = 3/2;
        res43(i) = 4/3;
        res54(i) = 5/4;
        res65(i) = 6/5;
        
        res31(i) = 3/1;
        res53(i) = 5/3;
        res75(i) = 7/5;
        
        res41(i) = 4/1;
        res52(i) = 5/2;
        res74(i) = 7/4;
        
        res12(i) = 1/2;
        res23(i) = 2/3;
        res34(i) = 3/4;
        res45(i) = 4/5;
        res56(i) = 5/6;
        
        res13(i) = 1/3;
        res35(i) = 3/5;
        res57(i) = 5/7;
        
        res14(i) = 1/4;
        res25(i) = 2/5;
        res47(i) = 4/7;
    end
    
    //Create figure
    
    fres = scf();
    fres.figure_name = "Resonances";
    fres.background = color(255,255,255);
    fres.figure_size = [1400,500];
    
    worldName = [string("Miranda");
             string("Ariel");
             string("Umbriel");
             string("Titania");
             string("Oberon")];

    for (m = 1:1:nmoons)
        subplot(1,nmoons,m);
        title(worldName(m),'fontsize',4);
        xlabel("Time (Gyr)",'fontsize',4);
        ylabel("Mean motion ratio",'fontsize',4);

       time5000 = 0:0.001:5;
        plot2d(time5000,res21);
        plot2d(time5000,res32);
        plot2d(time5000,res43);
        plot2d(time5000,res54);
        plot2d(time5000,res65);
        
        plot2d(time5000,res31,color("gray80"));
        plot2d(time5000,res53,color("gray80"));
        plot2d(time5000,res75,color("gray80"));
        
        plot2d(time5000,res41,color("gray90"));
        plot2d(time5000,res52,color("gray90"));
        plot2d(time5000,res74,color("gray90"));
        
        plot2d(time5000,res12);
        plot2d(time5000,res23);
        plot2d(time5000,res34);
        plot2d(time5000,res45);
        plot2d(time5000,res56);
        
        plot2d(time5000,res13,color("gray80"));
        plot2d(time5000,res35,color("gray80"));
        plot2d(time5000,res57,color("gray80"));
        
        plot2d(time5000,res14,color("gray90"));
        plot2d(time5000,res25,color("gray90"));
        plot2d(time5000,res47,color("gray90"));

        for (n = 1:1:nmoons)
            if (n <> m)
                aRatio = zeros(NT)
                for (t = 1:1:NT)
                    aRatio(t) = MMratio(m,n,t);
                end
                if (n == 1) then
                    plot2d(time(max(tspawn(m),tspawn(n)):NT),aRatio(max(tspawn(m),tspawn(n)):NT),color("deepskyblue"));
                    set_line_style(2.8);
                end
                if (n == 2) then
                    plot2d(time(max(tspawn(m),tspawn(n)):NT),aRatio(max(tspawn(m),tspawn(n)):NT),color("mediumseagreen"));
                    set_line_style(2.8);
                end
                if (n == 3) then
                    plot2d(time(max(tspawn(m),tspawn(n)):NT),aRatio(max(tspawn(m),tspawn(n)):NT),color("dark orange"));
                    set_line_style(2.8);
                end
                if (n == 4) then
                    plot2d(time(max(tspawn(m),tspawn(n)):NT),aRatio(max(tspawn(m),tspawn(n)):NT),color("darkviolet"));
                    set_line_style(2.8);
                end
                if (n == 5) then
                    plot2d(time(max(tspawn(m),tspawn(n)):NT),aRatio(max(tspawn(m),tspawn(n)):NT),color("deeppink"));
                    set_line_style(2.8);
                end
            end
        end
        a=get("current_axes") //get the handle of the newly created axes
        a.font_size = fontsize;
        a.axes_visible="on"; // makes the axes visible
        a.log_flags = "nn" ; // lin-lin
        a.data_bounds = [0,0;5,max(MMratio(m,:,:))];
        a.sub_ticks = [8 4];
    end
    
    //Export to PNG file
    xs2png(fres,"Res.png");
