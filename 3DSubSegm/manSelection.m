% Authors: Michał Ziemczonok, Martyna Mazur
% Contact: michal.ziemczonok.dokt@pw.edu.pl, martyna.mazur.dokt@pw.edu.pl
% Affiliation: Warsaw University of Technology, Institute of Micromechanics and Photonics

function f = manSelection(REC,seedAut,REC2,scale)
%input parsing - allow calls with any combination of arguments
if nargin == 1
    REC2 = [];
    scale = [];
    seedAut = [];
elseif nargin == 2
    REC2 = [];
    scale = [];
elseif nargin == 3
    if length(REC2(:)) == length(REC(:))
        scale = [];
    else
        scale = REC2;
        REC2 = [];
    end
end

%if variable is complex -> show log(abs())
if ~isreal(REC)
    disp('complex array - displaying log(0.1+abs(REC))');
    REC = log(0.1+abs(REC));
    REC2 = log(0.1+abs(REC2));
    if isempty(scale) %its probably bad anyway
        scale = [-2 12];
    end
end

%adjust the scale automatically
if isempty(scale) 
    scale = [REC(abs(REC(:)-median(REC(:))) > 0.005); REC2(abs(REC2(:)-median(REC2(:))) > 0.005)]; %account for data in both RECs
    scale = sort(scale); 
    scale = [scale(round(0.001*end))-0.008 scale(round(0.999*end))+0.008];
end

%initialize figure
f = figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [-0.0042, 0.0296, 1.0083, 0.9778]);
z = round(size(REC,3)/2); %default Z crossection
if isempty(REC2)
    imagesc(REC(:,:,z),scale); axis image;
    set(gca,'YDir', 'reverse'); %'normal' - reconstructions are upside down
    pos = get(gca,'Position');
    cb = colorbar;
    set(gca,'Position',pos);
else
    ax = subplot(1,2,1,'align');
    imagesc(REC(:,:,z),scale); axis image;
    ax.YDir = 'reverse';
    
    ax = subplot(1,2,2,'align');
    imagesc(REC2(:,:,z),scale); axis image;
    ax.YDir = 'reverse';
    
    pos = get(gca,'Position');
    cb = colorbar;
    set(gca,'Position',pos);
end

PushButton = uicontrol(gcf,'Style', 'push', 'String', 'Save Seed Points','Position', [200 100 200 60],'CallBack', @PushB,'FontSize',12);
% PushButton = uicontrol(gcf,'Style', 'push', 'String', 'Select Seed Points','Position', [200 200 200 60],'CallBack', @PushB,'FontSize',12);
             
seed=zeros(1,3);
diameter=0;
udata.start = [];
udata.axisNames = ['Y'; 'X'; 'Z'];
udata.changeNames = 1;
udata.sensitivity = 1; %how fast the planes are scrolled
udata.ZOffset = round(size(REC,3)/2);
udata.ZOffset2 = 0;
udata.REC = REC; clear REC;
udata.REC2 = REC2; clear REC2;
set(f,'UserData',udata);
set(f,'WindowButtonDown',@mouseDown);
set(f, 'WindowKeyPressFcn', @keyPressCallback);
set(f,'WindowButtonUp',@mouseUp);
set(f,'WindowButtonMotionFcn', @mouseMove);
set(f,'WindowScrollWheelFcn',@Scroll);

function mouseDown(object, eventdata) 
    f = gcf;
    
    if strcmp(f.SelectionType, 'normal') %LMB - change the plane
        f.UserData.start = get(f, 'CurrentPoint');
               
    elseif strcmp(f.SelectionType, 'alt') %RMB - "rotate" the reconstruction
        f.UserData.changeNames = 1; %improves the performance by a LOT
        %make sure the axis directions make intuitive sense
        if f.UserData.axisNames(3) == 'Z' %z = -z
            f.UserData.REC = permute(f.UserData.REC,[3 1 2]);
            f.UserData.REC2 = permute(f.UserData.REC2,[3 1 2]);
            f.UserData.axisNames = circshift(f.UserData.axisNames,1);
        elseif f.UserData.axisNames(3) == 'X' %swap Z and Y axis
            f.UserData.REC = permute(f.UserData.REC,[1 3 2]);
            f.UserData.REC2 = permute(f.UserData.REC2,[1 3 2]);
            f.UserData.axisNames = circshift(f.UserData.axisNames,1);
            f.UserData.axisNames(1:2) = circshift(f.UserData.axisNames(1:2),1);
        elseif f.UserData.axisNames(3) == 'Y' %revert z and x case changes
            f.UserData.REC = permute(f.UserData.REC,[3 2 1]);
            f.UserData.REC2 = permute(f.UserData.REC2,[3 2 1]);
            f.UserData.axisNames(1:2) = circshift(f.UserData.axisNames(1:2),1);
            f.UserData.axisNames = circshift(f.UserData.axisNames,1);
        end
    end
    f.UserData.sensitivity = (20+size(f.UserData.REC,3))/1;

end

function keyPressCallback(source,eventdata)
      f=gcf;
      ax = f.Children(3);
      cp = ax.CurrentPoint;

      keyPressed = eventdata.Key;
      if strcmp(keyPressed,'s')
         
          if f.UserData.axisNames(3) == 'Z'
           mousePos = ax.CurrentPoint;
           shift = cp(2) - mousePos(2);
           Zindex = round(f.UserData.ZOffset + shift*f.UserData.sensitivity);
               if Zindex == 0
                   Zindex = round(size(udata.REC,3)/2);
               end
             if ~isempty(seedAut)
                 hold on; h1 = plot(seedAut(:,1),seedAut(:,2),...
                'LineStyle',"none", 'Marker',"o", 'MarkerEdgeColor',"r",...
                'MarkerFaceColor',"r", 'MarkerSize', 3);
             end
             
            condition = true;
            h_p=[];
            
            while condition
                
                [Xseed,Yseed]=ginput(1);       
                new_seedPoints(:,1,:)=Xseed;
                new_seedPoints(:,2,:)=Yseed;
                new_seedPoints(:,3,:)=Zindex;

                seed_temp={seed, new_seedPoints};
                seed=vertcat(seed_temp{:});

                hold on; h2=plot(new_seedPoints(:,1,:),new_seedPoints(:,2,:),...
                    'LineStyle',"none", 'Marker',"o",'MarkerEdgeColor',"k",...
                    'MarkerFaceColor',"y", 'MarkerSize', 4);
                    h_p(end+1)=h2;
                    
                    if strcmp(f.SelectionType,'alt')
                        condition=false;
                        delete(h1);
                        delete(h_p);
%                         Z_seed=input('Z = ');
%                         seed(ismember(seed(:,3),0),3)=Z_seed;
                        seed(end,:)=[];
                    end
            
            end
            
             hold on; h3=plot(seed(:,1),seed(:,2),...
                'LineStyle',"none", 'Marker',"o",'MarkerEdgeColor',"k",...
                'MarkerFaceColor',"y", 'MarkerSize', 4);

            clear new_seedPoints;
            w = waitforbuttonpress;
%             delete(h1);
            delete(h3);
            
          end 
          
      elseif strcmp(keyPressed,'r')
            con = 0;
            XY = zeros(1,2);
            hall = [];
            while con ~= 2          
            [X_radius,Y_radius]=ginput(1);
            XY_temp={XY, [X_radius, Y_radius]};
            XY=vertcat(XY_temp{:});
           
            con = con + size(X_radius,1);
            hold on; hl2=plot(X_radius(:,1),Y_radius(:,1),...
                    'LineStyle',"none", 'Marker',"o",'MarkerEdgeColor',"k",...
                    'MarkerFaceColor',"m", 'MarkerSize', 4);
            hall(end+1) = hl2;
            
            end
            XY(1,:) = [];
            diameter = sqrt(((XY(1,1)-XY(2,1))^2)+((XY(1,2)-XY(2,2))^2));
            hold on;
            hl = line('XData',XY(:,1),'YData',XY(:,2),...
             'Marker','o','color','m');
            w = waitforbuttonpress;
            delete(hl);
            delete(hall)
            
            
      elseif strcmp(keyPressed,'escape')
            close(f);
          
      end
  end

function mouseUp(object, eventdata)
    f = gcf;
    f.UserData.start = [];
    %save the offset to start scrolling from the plane when the mouse was released
    try
        f.UserData.ZOffset = f.UserData.ZOffset2;
    catch error
    end
end

function mouseMove (object, eventdata)
    f = gcf;
    if ~isempty(f.UserData.start)
        mousePos = get(f, 'CurrentPoint');
        shift = f.UserData.start(2) - mousePos(2); %diff [ranging from 0 to 1] in figure coordinates
        
        %start scrolling from the middle plane and clear the offset
        if f.UserData.changeNames == 1
            f.UserData.ZOffset = round(size(f.UserData.REC,3)/2);
            f.UserData.ZOffset2 = 0;
        end

        %scale to Z values
        Zindex = round(f.UserData.ZOffset + shift*f.UserData.sensitivity);
        %update the last Z index
        f.UserData.ZOffset2 = Zindex;
        %min/max boundary check
        if Zindex >= size(f.UserData.REC,3); Zindex = size(f.UserData.REC,3); end
        if Zindex < 1; Zindex = 1; end
        
        %update figure, single input case
        if isempty(f.UserData.REC2)
            ax = f.Children(3);
            target = get(ax,'Children');
            set(target,'CData',f.UserData.REC(:,:,Zindex));
            set(get(ax, 'title'),'string',[f.UserData.axisNames(3) ' = ' num2str(Zindex)]);
            if f.UserData.changeNames == 1
                ax.YLim = [0.5 size(f.UserData.REC,1)];
                ax.XLim = [0.5 size(f.UserData.REC,2)];
                set(get(ax, 'ylabel'),'string',f.UserData.axisNames(1));
                set(get(ax, 'xlabel'),'string',f.UserData.axisNames(2));
                f.UserData.changeNames = 0;
            end

        else %comparison case
            axLeft = f.Children(4);
            target1 = get(axLeft,'Children');
            set(target1,'CData',f.UserData.REC(:,:,Zindex));
            set(get(axLeft, 'title'),'string',[f.UserData.axisNames(3) ' = ' num2str(Zindex)]);
            if f.UserData.changeNames == 1
                axLeft.YLim = [0.5 size(f.UserData.REC,1)];
                axLeft.XLim = [0.5 size(f.UserData.REC,2)];
                set(get(axLeft, 'ylabel'),'string',f.UserData.axisNames(1));
                set(get(axLeft, 'xlabel'),'string',f.UserData.axisNames(2));
            end

            axRight = f.Children(3);
            target2 = get(axRight,'Children');
            set(target2,'CData',f.UserData.REC2(:,:,Zindex));
            set(get(axRight, 'title'),'string',[f.UserData.axisNames(3) ' = ' num2str(Zindex)]);
            if f.UserData.changeNames == 1
                axRight.YLim = [0.5 size(f.UserData.REC,1)];
                axRight.XLim = [0.5 size(f.UserData.REC,2)];
                set(get(axRight, 'ylabel'),'string',f.UserData.axisNames(1));
                set(get(axRight, 'xlabel'),'string',f.UserData.axisNames(2));
                f.UserData.changeNames = 0;
            end
        end
    end
    pause(0.01);
end

%mouse scroll adjusts the zoom of the slice
function Scroll(object, eventdata)
    f = gcf;

    sensitivity = 12;
    
    y_limits = size(f.UserData.REC,1);
    x_limits = size(f.UserData.REC,2);
    mousePos = get(f, 'CurrentPoint'); %cursor location
    y_weight = (0.5-mousePos(2));%speed up the zoom in the direction which is closer to the edge of the figure
    x_weight = (0.5-mousePos(1));
    rectangularSensitivity = y_limits/x_limits;
    ScrollCount = eventdata.VerticalScrollCount * sensitivity;

    ax = f.Children(3);
    try
    ax.YLim = [ax.YLim(1) - (0.5+y_weight)*ScrollCount*rectangularSensitivity ax.YLim(2) + (0.5-y_weight)*ScrollCount*rectangularSensitivity];
    ax.XLim = [ax.XLim(1) - (0.5-x_weight)*ScrollCount ax.XLim(2) + (0.5+x_weight)*ScrollCount];
    catch
    end
    %check the limits
    if ax.YLim(1) < 1; ax.YLim(1) = 0.5; end
    if ax.YLim(2) > y_limits; ax.YLim(2) = y_limits+0.5; end
    if ax.XLim(1) < 1; ax.XLim(1) = 0.5; end
    if ax.XLim(2) > x_limits; ax.XLim(2) = x_limits+0.5; end

    %case with 2 images
    if ~isempty(f.UserData.REC2)
        ax = f.Children(4);    
        try
        ax.YLim = [ax.YLim(1) - (0.5+y_weight)*ScrollCount*rectangularSensitivity ax.YLim(2) + (0.5-y_weight)*ScrollCount*rectangularSensitivity];
        ax.XLim = [ax.XLim(1) - (0.5-x_weight)*ScrollCount ax.XLim(2) + (0.5+x_weight)*ScrollCount];
        catch
        end
        if ax.YLim(1) < 1; ax.YLim(1) = 0.5; end
        if ax.YLim(2) > y_limits; ax.YLim(2) = y_limits+0.5; end
        if ax.XLim(1) < 1; ax.XLim(1) = 0.5; end
        if ax.XLim(2) > x_limits; ax.XLim(2) = x_limits+0.5; end
    end
end
      

 function PushB(object,eventdata)
     
        if sum(seed(1,:,:))==0
            seed(1,:,:)=[];
        end
        temp = seed;
        seed = {};
        seed{1} = temp;
        seed{2} = diameter;
        save(['.\seed.mat'],'seed','-v7.3')
        
        t=text((size(udata.REC,1))/2,(size(udata.REC,2))/2,'Seed points were saved!','horizontal','center','color','r','FontSize',15);
        pause(1);
        delete(t);
%         pause(0.5);
%         close(f)

 end
end