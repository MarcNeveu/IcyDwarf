//    clear;          // Clear variables
    xdel(winsid()); // Clear figures

    nmoons = 5;        //Number of moons
    NT = 149-10+1;    //Number of timesteps (incl. zeroth)
    timestep = 0.01;  //Timestep of the sim, specified in IcyDwarfInput.txt
    time = zeros(NT);
    
    tspawn = zeros(nmoons); // Time of formation (Myr)
    tspawn(1) = 424-10+1; // Mimas
    tspawn(2) = 354-10+1;  // Enceladus
    tspawn(3) = 150-10+1;  // Tethys
    tspawn(4) = 10-10+1;  // Dione
    tspawn(5) = 10-10+1;  // Rhea
    
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
    Tmin = 3000;
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
                 [num_read, Radius(m,t,r), Temp(m,t,r), Mrock(m,t,r), Mh2os(m,t,r), Madhs(m,t,r), Mh2ol(m,t,r), Mnh3l(m,t,r), Nu(m,t,r), Famor(m,t,r), Kappa(m,t,r), Xhydr(m,t,r), Pore(m,t,r)] = mfscanf(fid(m), "%f %f %f %f %f %f %f %f %f %f %f %f");
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
    nz = [50 100 150 200 250 300 400 500 600 700 800 900 1000 1100 1200 1300];

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
            title("Mimas",'color',color("deepskyblue"),'fontsize',6);
            ylabel("Radius (km)",'fontsize',6);
        end
        if (m==2)
             title("Enceladus",'color',color("spring green"),'fontsize',6);   
        end
        if (m==3)
             title("Tethys",'color',color("dark orange"),'fontsize',6);   
        end
        if (m==4)
             title("Dione",'color',color("darkviolet"),'fontsize',6);   
        end
        if (m==5)       
             title("Rhea",'color',color("deeppink"),'fontsize',6); 
        end
        xlabel("Time (Gyr)",'fontsize',6);
        
        Sgrayplot(time(tspawn(m):NT),Radiustemp,Temptemp(tspawn(m):NT,:),rect=[0,0,time(50),Radiustemp(NR)],zminmax=[0 2000]);
        
        xset('font size',6);
        //xset("fpf","%.0f");
        xset("fpf"," ");
        //xset('thickness',1.5);
        contour2d(time(tspawn(m):NT),Radiustemp,Temptemp(tspawn(m):NT,:),nz,rect=[0,0,time(50),Radiustemp(NR)]);
        
        a=get("current_axes") //get the handle of the newly created axes
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
                if (Mrock(m,t,r) <= 0 & Mrock(m,t,r-1) > 0)
                    Core(m,t) = Radius(m,t,r)
                end
                if (Mrock(m,t,r) > 0 & Mrock(m,t,r-1) <= 0 & Mh2ol(m,t,r-1) <= 0)
                    Ice(m,t) = Radius(m,t,r)
                end
                if (Mrock(m,t,r) <= 0 & Mh2ol(m,t,r) > 0)
                    Ocean(m,t) = Radius(m,t,r)
                end
                if (Mrock(m,t,r) > 0 & Mh2ol(m,t,r) > 0)
                    Liqcore(m,t) = Radius(m,t,r)
                end
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
            title("Mimas",'color',color("deepskyblue"),'fontsize',6);
            ylabel("Radius (km)",'fontsize',6);
        end
        if (m==2)
             title("Enceladus",'color',color("spring green"),'fontsize',6);   
        end
        if (m==3)
             title("Tethys",'color',color("dark orange"),'fontsize',6);   
        end
        if (m==4)
             title("Dione",'color',color("darkviolet"),'fontsize',6);   
        end
        if (m==5)       
             title("Rhea",'color',color("deeppink"),'fontsize',6); 
        end
        xlabel("Time (Gyr)",'fontsize',6);
        
        xfpoly(timepoly(m,tspawn(m):$),Crust(m,tspawn(m):$)); // ' = transpose
        h = gce();
        h.line_mode = "off";
        h.background = color("gray50"); // Adjust color
        
        if (max(Ice(m,tspawn(m):$)) > 0) then
            xfpoly(timepoly(m,tspawn(m):$),Ice(m,tspawn(m):$)); // ' = transpose
            h = gce();
            h.line_mode = "off";
            h.background = color("lightblue2"); // Adjust color
        end
        
        if (max(Ocean(m,tspawn(m):$)) > 0) then
            xfpoly(timepoly(m,tspawn(m):$),Ocean(m,tspawn(m):$)); // ' = transpose
            h = gce();
            h.line_mode = "off";
            h.background = color("royalblue1"); // Adjust color
        end
        
        if (max(Core(m,tspawn(m):$)) > 0) then
            xfpoly(timepoly(m,tspawn(m):$),Core(m,tspawn(m):$)); // ' = transpose
            h = gce();
            h.line_mode = "off";
            h.background = color("chocolate"); // Adjust color
        end
        
        if (max(Liqcore(m,tspawn(m):$)) > 0) then
            xfpoly(timepoly(m,tspawn(m):$),Liqcore(m,tspawn(m):$)); // ' = transpose
            h = gce();
            h.line_mode = "off";
            h.background = color("light sea green"); // Adjust color
        end
        
        xset('font size',6);
        xset("fpf"," ");
        
        a=get("current_axes") //get the handle of the newly created axes
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
    
//    NT = NT+10-1
//    NT = NT*10;
//    NT = NT-100+1;
//    for (m = 1:1:nmoons)
//        tspawn(m) = tspawn(m)+10-1;
//        tspawn(m) = tspawn(m)*10;
//        tspawn(m) = tspawn(m)-100+1;
//        if (tspawn(i) > NT) 
//            tspawn(i) = NT-5;
//        end
//    end
//    timestep = timestep/10;
//    
//    time = zeros(NT);
//    for (m=1:1:nmoons)
//        for (t = 1:1:NT) 
//            if (t>1) 
//                time(t) = time(t-1) + timestep;
//            end
//        end
//    end
//    
//    TITLEOrb = string("");
//    infile(1) = '0Orbit.txt';
//    infile(2) = '1Orbit.txt';
//    infile(3) = '2Orbit.txt';
//    infile(4) = '3Orbit.txt';
//    infile(5) = '4Orbit.txt';
//    
//    timeorb = zeros(NT);
//    aorb = zeros(nmoons,NT);
//    aorbtemp = zeros(NT);
//    aosc = zeros(nmoons,NT);
//    eorb = zeros(nmoons,NT);
//    eorbtemp = zeros(NT);
//    horb = zeros(nmoons,NT);
//    korb = zeros(nmoons,NT);
//    resangle = zeros(nmoons,NT);
//    Wtide = zeros(nmoons,NT);
//    Wtidetemp = zeros(NT);
//    k2Q = zeros(nmoons,NT);
//    k2Qtemp = zeros(NT);
//        
//    mooncolor = zeros(nmoons);
//    
//    mooncolor(1) = color("deepskyblue");
//    mooncolor(2) = color("spring green");
//    mooncolor(3) = color("dark orange");
//    mooncolor(4) = color("darkviolet");
//    mooncolor(5) = color("deeppink");
//    
//    moonsize(1) = 198;
//    moonsize(2) = 252;
//    moonsize(3) = 530;
//    moonsize(4) = 561;
//    moonsize(5) = 763;
//
//    //Read the orbital output files
//    
//    fidorb = zeros(nmoons);
//    
//    for (m = 1:1:nmoons)
//        fidorb(m) = mopen(infile(m), "r");
//        if (fidorb(m) == -1) then
//            error("Cannot open input file");
//        end
//        
//        for (t = 1:1:NT)
//            [num_read, timeorb(t), aorb(m,t), aosc(m,t), eorb(m,t), horb(m,t), korb(m,t), resangle(m,t), Wtide(m,t), k2Q(m,t)] = mfscanf(fidorb(m), "%f %f %f %f %f %f %f %f %f");
//        end   
//        
//        mclose(fidorb(m)); 
//    end
//
//    //Create figure
//    
//    forb = scf();
//    forb.figure_name = "Orbits";
//    forb.background = color(255,255,255);
//    forb.figure_size = [1200,600];
//
//    // Semi-major axes
//    subplot(1,4,1);
//    xset('font size',6);
//    xset("fpf"," ");
//    xset('thickness',2.8);
//    title("Semi-major axes",'fontsize',6);
//    xlabel("Time (Gyr)",'fontsize',6);
//    ylabel("Semi-major axis (km)",'fontsize',6);
//    for (m = 1:1:nmoons)
//        for (t = 1:1:NT)
//            aorbtemp(t) = aorb(m,t);
//        end
//        plot2d(time(tspawn(m):NT),aorbtemp(tspawn(m):NT),mooncolor(m));
//    end
//    a=get("current_axes") //get the handle of the newly created axes
//    a.axes_visible="on"; // makes the axes visible
//    a.log_flags = "nn" ; // lin-lin
//    a.data_bounds(:,1) = [0;5] ; // set positive bounds for X axis
//    a.sub_ticks = [8 4];
//    
//    // Overlay current semi-major axes as filled circles
//    currentaorb = zeros(nmoons);
//    currentaorb(1) = 185539;
//    currentaorb(2) = 237948;
//    currentaorb(3) = 294660;
//    currentaorb(4) = 377400;
//    currentaorb(5) = 527040;
//    
//    currenttime = zeros(nmoons);
//    for (m = 1:1:nmoons)
//        currenttime(m) = 4.56;
//    end
//    
//    if (nmoons==5)
//        scatter(currenttime,currentaorb,moonsize,[color("deepskyblue")+1;color("spring green")+1;color("dark orange")+1;4;color("deeppink")],"fill",".");
//    end
//    if (nmoons <> 5) then
//        mprintf("Not plotting the current orbital positions because nmoons < 5\n");
//    end
//
//    // Eccentricities
//    subplot(1,4,2); 
//    xset('font size',6);
//    xset('thickness',2.8);
//    title("Eccentricities",'fontsize',6);
//    xlabel("Time (Gyr)",'fontsize',6);
//    ylabel("Eccentricity",'fontsize',6);
//    for (m = 1:1:nmoons)
//        for (t = 1:1:NT)
//            eorbtemp(t) = eorb(m,t);
//        end
//        plot2d(time(tspawn(m):NT),eorbtemp(tspawn(m):NT),mooncolor(m));//,rect=[0,1e-7,time(50),1]);
//    end
//    a=get("current_axes") //get the handle of the newly created axes
//    a.axes_visible="on"; // makes the axes visible
//    a.log_flags = "nl" ; // lin-log
//    a.data_bounds = [0,1e-7;5,1];
//    a.sub_ticks = [8 4];
//    
//    // Overlay current eccentricities as filled circles
//    currenteorb = zeros(nmoons);
//    currenteorb(1) = 0.0196;
//    currenteorb(2) = 0.0047;
//    currenteorb(3) = 0.0001;
//    currenteorb(4) = 0.0022;
//    currenteorb(5) = 0.0010;
//    
//    if (nmoons==5)
//        scatter(currenttime,currenteorb,moonsize,[color("deepskyblue")+1;color("spring green")+1;color("dark orange")+1;4;color("deeppink")],"fill",".");
//    end
//    if (nmoons <> 5) then
//        mprintf("Not plotting the current orbital positions because nmoons < 5\n");
//    end
//    
//    // Tidal dissipation
//    subplot(1,4,3); 
//    xset('font size',6);
//    xset('thickness',2.8);
//    title("Total tidal dissipation",'fontsize',6);
//    xlabel("Time (Gyr)",'fontsize',6);
//    ylabel("Tidal dissipation (W)",'fontsize',6);
//    for (m = 1:1:nmoons)
//        for (t = 1:1:NT)
//            Wtidetemp(t) = Wtide(m,t);
//        end
//        plot2d(time(tspawn(m):NT),Wtidetemp(tspawn(m):NT),mooncolor(m));//,rect=[0,1e-7,time(50),1]);
//    end
//    a=get("current_axes") //get the handle of the newly created axes
//    a.axes_visible="on"; // makes the axes visible
//    a.data_bounds = [0,1e-1;5,1e12];
//    a.log_flags = "nl" ; // lin-log
//    a.sub_ticks = [8 4];
//    
//    // Equivalent k2/Q
//    subplot(1,4,4); 
//    xset('font size',6);
//    xset('thickness',2.8);
//    title("k2/Q",'fontsize',6);
//    xlabel("Time (Gyr)",'fontsize',6);
//    ylabel("k2/Q",'fontsize',6);
//    for (m = 1:1:nmoons)
//        for (t = 1:1:NT)
//            k2Qtemp(t) = k2Q(m,t);
//        end
//        plot2d(time(tspawn(m):NT),k2Qtemp(tspawn(m):NT),mooncolor(m));//,rect=[0,1e-7,time(50),1]);
//    end
//    a=get("current_axes") //get the handle of the newly created axes
//    a.axes_visible="on"; // makes the axes visible
//    a.data_bounds = [0,1e-12;5,1e-2];
//    a.log_flags = "nl" ; // lin-log
//    a.sub_ticks = [8 4];
//    
//    //Export to PNG file
//    xs2png(forb,"Orb.png");
//    
//    // -----------------------------------
//    // RESONANCES
//    // -----------------------------------
//
//    nME = zeros(NT);
//    nMT = zeros(NT);
//    nMD = zeros(NT);
//    nMR = zeros(NT);
//    nET = zeros(NT);
//    nED = zeros(NT);
//    nER = zeros(NT);
//    nTD = zeros(NT);
//    nTR = zeros(NT);
//    nDR = zeros(NT);
//    
//    for (i = 1:1:NT)
//        nME(i) = (aorb(1,i)/aorb(2,i))^(-1.5);
//        nMT(i) = (aorb(1,i)/aorb(3,i))^(-1.5);
//        nMD(i) = (aorb(1,i)/aorb(4,i))^(-1.5);
//        nMR(i) = (aorb(1,i)/aorb(5,i))^(-1.5);
//        nET(i) = (aorb(2,i)/aorb(3,i))^(-1.5);
//        nED(i) = (aorb(2,i)/aorb(4,i))^(-1.5);
//        nER(i) = (aorb(2,i)/aorb(5,i))^(-1.5);
//        nTD(i) = (aorb(3,i)/aorb(4,i))^(-1.5);
//        nTR(i) = (aorb(3,i)/aorb(5,i))^(-1.5);
//        nDR(i) = (aorb(4,i)/aorb(5,i))^(-1.5);
//    end
//    
//    // Overlay resonances
//    res21 = zeros(5001);
//    res32 = zeros(5001);
//    res43 = zeros(5001);
//    res54 = zeros(5001);
//    res65 = zeros(5001);
//    
//    res31 = zeros(5001);
//    res53 = zeros(5001);
//    res75 = zeros(5001);
//    
//    res41 = zeros(5001);
//    res52 = zeros(5001);
//    res74 = zeros(5001);
//    
//    for (i = 1:1:5001)
//        res21(i) = 2/1;
//        res32(i) = 3/2;
//        res43(i) = 4/3;
//        res54(i) = 5/4;
//        res65(i) = 6/5;
//        
//        res31(i) = 3/1;
//        res53(i) = 5/3;
//        res75(i) = 7/5;
//        
//        res41(i) = 4/1;
//        res52(i) = 5/2;
//        res74(i) = 7/4;
//    end
//    
//    //Create figure
//    
//    fres = scf();
//    fres.figure_name = "Orbits";
//    fres.background = color(255,255,255);
//    fres.figure_size = [1200,600];
//
//    // Mimas
//    subplot(1,4,1);
//    xset('font size',6);
//    xset("fpf"," ");
//    xset('thickness',2.8);
//    title("Mimas",'color',mooncolor(1),'fontsize',6);
//    xlabel("Time (Gyr)",'fontsize',6);
//    ylabel("Mean motion ratio",'fontsize',6);
//
//    plot2d(time(tspawn(1):NT),nME(tspawn(1):NT),mooncolor(2));
//    plot2d(time(tspawn(1):NT),nMT(tspawn(1):NT),mooncolor(3));
//    plot2d(time(tspawn(1):NT),nMD(tspawn(1):NT),mooncolor(4));
//    plot2d(time(tspawn(1):NT),nMR(tspawn(1):NT),mooncolor(5));
//    a=get("current_axes") //get the handle of the newly created axes
//    a.axes_visible="on"; // makes the axes visible
//    a.log_flags = "nn" ; // lin-lin
//    a.data_bounds(:,1) = [0;5] ; // set positive bounds for X axis
//    a.sub_ticks = [8 4];
//    
//    time5000 = 0:0.001:5;
//    plot2d(time5000,res21);
//    plot2d(time5000,res32);
//    plot2d(time5000,res43);
//    plot2d(time5000,res54);
//    plot2d(time5000,res65);
//        
//    plot2d(time5000,res31,color("gray80"));
//    plot2d(time5000,res53,color("gray80"));
//    plot2d(time5000,res75,color("gray80"));
//    
//    plot2d(time5000,res41,color("gray90"));
//    plot2d(time5000,res52,color("gray90"));
//    plot2d(time5000,res74,color("gray90"));
//    
//    // Enceladus
//    subplot(1,4,2);
//    xset('font size',6);
//    xset("fpf"," ");
//    xset('thickness',2.8);
//    title("Enceladus",'color',mooncolor(2),'fontsize',6);
//    xlabel("Time (Gyr)",'fontsize',6);
//    ylabel("Mean motion ratio",'fontsize',6);
//
//    plot2d(time(tspawn(2):NT),nET(tspawn(2):NT),mooncolor(3));
//    plot2d(time(tspawn(2):NT),nED(tspawn(2):NT),mooncolor(4));
//    plot2d(time(tspawn(2):NT),nER(tspawn(2):NT),mooncolor(5));
//    a=get("current_axes") //get the handle of the newly created axes
//    a.axes_visible="on"; // makes the axes visible
//    a.log_flags = "nn" ; // lin-lin
//    a.data_bounds(:,1) = [0;5] ; // set positive bounds for X axis
//    a.sub_ticks = [8 4];
//    
//    plot2d(time5000,res21);
//    plot2d(time5000,res32);
//    plot2d(time5000,res43);
//    plot2d(time5000,res54);
//    plot2d(time5000,res65);
//    
//    plot2d(time5000,res31,color("gray80"));
//    plot2d(time5000,res53,color("gray80"));
//    plot2d(time5000,res75,color("gray80"));
//    
//    plot2d(time5000,res41,color("gray90"));
//    plot2d(time5000,res52,color("gray90"));
//    plot2d(time5000,res74,color("gray90"));
//    
//    // Tethys
//    subplot(1,4,3);
//    xset('font size',6);
//    xset("fpf"," ");
//    xset('thickness',2.8);
//    title("Tethys",'color',mooncolor(3),'fontsize',6);
//    xlabel("Time (Gyr)",'fontsize',6);
//    ylabel("Mean motion ratio",'fontsize',6);
//
//    plot2d(time(tspawn(3):NT),nTD(tspawn(3):NT),mooncolor(4));
//    plot2d(time(tspawn(3):NT),nTR(tspawn(3):NT),mooncolor(5));
//    a=get("current_axes") //get the handle of the newly created axes
//    a.axes_visible="on"; // makes the axes visible
//    a.log_flags = "nn" ; // lin-lin
//    a.data_bounds(:,1) = [0;5] ; // set positive bounds for X axis
//    a.sub_ticks = [8 4];
//    
//    plot2d(time5000,res21);
//    plot2d(time5000,res32);
//    plot2d(time5000,res43);
//    plot2d(time5000,res54);
//    plot2d(time5000,res65);
//    
//    plot2d(time5000,res31,color("gray80"));
//    plot2d(time5000,res53,color("gray80"));
//    plot2d(time5000,res75,color("gray80"));
//    
//    plot2d(time5000,res41,color("gray90"));
//    plot2d(time5000,res52,color("gray90"));
//    plot2d(time5000,res74,color("gray90"));
//    
//    // Dione
//    subplot(1,4,4);
//    xset('font size',6);
//    xset("fpf"," ");
//    xset('thickness',2.8);
//    title("Dione",'color',mooncolor(4),'fontsize',6);
//    xlabel("Time (Gyr)",'fontsize',6);
//    ylabel("Mean motion ratio",'fontsize',6);
//
//    plot2d(time(tspawn(4):NT),nDR(tspawn(4):NT),mooncolor(5));
//    a=get("current_axes") //get the handle of the newly created axes
//    a.axes_visible="on"; // makes the axes visible
//    a.log_flags = "nn" ; // lin-lin
//    a.data_bounds(:,1) = [0;5] ; // set positive bounds for X axis
//    a.sub_ticks = [8 4];
//    
//    plot2d(time5000,res21);
//    plot2d(time5000,res32);
//    plot2d(time5000,res43);
//    plot2d(time5000,res54);
//    plot2d(time5000,res65);
//        
//    plot2d(time5000,res31,color("gray80"));
//    plot2d(time5000,res53,color("gray80"));
//    plot2d(time5000,res75,color("gray80"));
//    
//    plot2d(time5000,res41,color("gray90"));
//    plot2d(time5000,res52,color("gray90"));
//    plot2d(time5000,res74,color("gray90"));
//    
//    //Export to PNG file
//    xs2png(fres,"Res.png");
