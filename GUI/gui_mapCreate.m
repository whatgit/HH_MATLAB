%{
mapCreate function use to create a map of warehouse environment
type :
    - 4tags : straight line dataset (first record)
    - 8tags_cw  : clockwise dataset (first record) 
    - 8tags_ccw : counter-clockwise dataset (first record)
    - start_tag_1   : start from tag no. 1 ; dataset 1,2,6 (second record)
    - start_tag_28  : start fron tag no. 28 ; dataset 8,10 (second record)
    - start_tag_21  : start fron tag no. 21 ; dataset 3 (second record)
    - start_tag_22  : start from tag no. 22 ; dataset 7 (second record)
    - start_tag_16  : start from tag np. 16 ; dataset 5 (second record)
%}

function gui_mapCreate(type)

    tag_scale = 15;  %15x tag size
    
    first_R_pillar_x = (2.94/2);
    first_R_pillar_y = 0;
    
    first_L_pillar_x = -(2.94/2);
    first_L_pillar_y = 0;
    
    gap_vertical = 2.70;    %vertical gap between pillars
    gap_horizontal = 1.77;  %gap in the middle between two row of pillars
    
    img_go_up = imread('tag.jpg'); 
    img_go_right = imrotate(img_go_up, 90);
    img_go_down = imrotate(img_go_up, 179);
    img_go_left = imrotate(img_go_up, 270);
    
    if strcmp(type,'4tags')
    
        %Plot pillars        
            %Left Pillars
            plot(-1.425,1,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-1.425 -1.425],[1 3.7],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-1.425 -1.425],[3.7 6.4],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-1.425 -1.425],[6.4 9.1],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-1.425 -1.425],[9.1 11.8],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
            %Right Pillars    
            plot(1.425,1,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([1.425 1.425],[1 3.7],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([1.425 1.425],[3.7 6.4],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([1.425 1.425],[6.4 9.1],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([1.425 1.425],[9.1 11.8],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);

    end
    
    if strcmp(type,'8tags_cw')
            
        x_fixer = 4.62;
        y_fixer = 4.0;
        
        %Plot tags
        %A4 Paper size 297 x 210 mm   
        %image([x vector],[y vector],image]
        image([-1.025+x_fixer -1.025+(0.0210*tag_scale)+x_fixer],[0.15+y_fixer 0.15+(0.0297*tag_scale)+y_fixer],img_go_up); %tags#1
        image([0.935+x_fixer 0.935+(0.0210*tag_scale)+x_fixer],[2.85+y_fixer 2.85+(0.0297*tag_scale)+y_fixer],img_go_up);   %tags#2
        image([0.265+x_fixer-(0.0297*tag_scale) 0.265+x_fixer],[4.25+y_fixer 4.25+(0.0210*tag_scale)+y_fixer],img_go_left); %tags#3
        image([-1.935+x_fixer-(0.0297*tag_scale) -1.935+x_fixer],[4.74+y_fixer 4.74+(0.0210*tag_scale)+y_fixer],img_go_left); %tags#4
        image([-4.505+x_fixer-(0.0210*tag_scale) -4.505+x_fixer],[4.74+y_fixer-(0.0297*tag_scale) 4.74+y_fixer],img_go_down); %tags#5
        image([-5.385+x_fixer-(0.0210*tag_scale) -5.385+x_fixer],[3.28+y_fixer-(0.0297*tag_scale) 3.28+y_fixer],img_go_down); %tags#6
        image([-3.835+x_fixer -3.835+(0.0210*tag_scale)+x_fixer],[0.58+y_fixer-(0.0297*tag_scale) 0.58+y_fixer],img_go_down); %tags#7
        %image([-4.765-(0.0210*tag_scale) -4.765],[-2.12-(0.0297*tag_scale) -2.12],img_go_down); %tags#8
        image([-4.600+x_fixer-(0.0210*tag_scale) -4.600+x_fixer],[-2.12+y_fixer-(0.0297*tag_scale) -2.12+y_fixer],img_go_down); %tags#8
    
        %Plot pillars        
        %Corridor A
            %Left Pillars
            plot([-6.045+x_fixer -6.045+x_fixer],[3.0+y_fixer 5.7+y_fixer] ,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-6.045+x_fixer -6.045+x_fixer],[0.3+y_fixer 3.0+y_fixer] ,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-6.045+x_fixer -6.045+x_fixer],[0.3+y_fixer -3.0+y_fixer] ,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
            %Right Pillars   
            
            plot([-1.425+x_fixer -3.195+x_fixer],[5.7+y_fixer 5.7+y_fixer],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            
            plot([-1.425+x_fixer -3.195+x_fixer],[3.0+y_fixer 3.0+y_fixer],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-3.195+x_fixer -3.195+x_fixer],[3.0+y_fixer 0.3+y_fixer],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-3.195+x_fixer -3.195+x_fixer],[0.3+y_fixer -3.0+y_fixer],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
        
        %Corridor B
            %Left Pillars
            plot(-1.425+x_fixer,0.3+y_fixer,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-1.425+x_fixer -1.425+x_fixer],[0.3+y_fixer 3.0+y_fixer],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot(-1.425+x_fixer,5.70+y_fixer,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
            %Right Pillars    
            plot(1.425+x_fixer,0.3+y_fixer,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([1.425+x_fixer 1.425+x_fixer],[0.3+y_fixer 3.0+y_fixer],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([1.425+x_fixer 1.425+x_fixer],[3.0+y_fixer 5.7+y_fixer],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
        
    end

    if strcmp(type,'8tags_ccw')
   
        %Plot tags
        %A4 Paper size 297 x 210 mm   
        %image([x vector],[y vector],image]
        image([-1.025 -1.025+(0.0210*tag_scale)],[0.15 0.15+(0.0297*tag_scale)],img_go_up); %tags#1
        image([0.935 0.935+(0.0210*tag_scale)],[2.85 2.85+(0.0297*tag_scale)],img_go_up);   %tags#2
        image([0.265-(0.0297*tag_scale) 0.265],[4.25 4.25+(0.0210*tag_scale)],img_go_left); %tags#3
        image([-1.935-(0.0297*tag_scale) -1.935],[4.74 4.74+(0.0210*tag_scale)],img_go_left); %tags#4
        image([-4.505-(0.0210*tag_scale) -4.505],[4.74-(0.0297*tag_scale) 4.74],img_go_down); %tags#5
        image([-5.385-(0.0210*tag_scale) -5.385],[3.28-(0.0297*tag_scale) 3.28],img_go_down); %tags#6
        image([-3.835 -3.835+(0.0210*tag_scale)],[0.58-(0.0297*tag_scale) 0.58],img_go_down); %tags#7
        %image([-4.765-(0.0210*tag_scale) -4.765],[-2.12-(0.0297*tag_scale) -2.12],img_go_down); %tags#8
        image([-4.600-(0.0210*tag_scale) -4.600],[-2.12-(0.0297*tag_scale) -2.12],img_go_down); %tags#8
    
        %Plot pillars        
        %Corridor A
            %Left Pillars
            plot([-6.045 -6.045],[3.0 5.7] ,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-6.045 -6.045],[0.3 3.0] ,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-6.045 -6.045],[0.3 -3.0] ,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
            %Right Pillars   
            
            plot([-1.425 -3.195],[5.7 5.7],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            
            plot([-1.425 -3.195],[3.0 3.0],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-3.195 -3.195],[3.0 0.3],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-3.195 -3.195],[0.3 -3.0],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
        
        %Corridor B
            %Left Pillars
            plot(-1.425,0.3,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([-1.425 -1.425],[0.3 3.0],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot(-1.425,5.70,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
            %Right Pillars    
            plot(1.425,0.3,'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([1.425 1.425],[0.3 3.0],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
            plot([1.425 1.425],[3.0 5.7],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
    end     
    
    if strcmp(type,'start_tag_1')    %Straight ccw start from big tag 
        
%         plot(first_R_pillar_x,first_R_pillar_y,'--rs','LineWidth',2,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','g',...
%                     'MarkerSize',10);
%        plot(first_L_pillar_x,first_L_pillar_y,'--rs','LineWidth',2,...
%             'MarkerEdgeColor','k',...
%             'MarkerFaceColor','g',...
%             'MarkerSize',10);

       %Corridor B 
       
       %axis([-7 2 -5 42]);
       
       for i=1:12
           
           plot([first_R_pillar_x,first_R_pillar_x],[first_R_pillar_y,first_R_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
           plot([first_L_pillar_x,first_L_pillar_x],[first_L_pillar_y,first_L_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
           
       end
       
       %Corridor A
       for i=1:12
           
           plot([first_L_pillar_x-gap_horizontal,first_L_pillar_x-gap_horizontal],[first_R_pillar_y,first_R_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
           plot([first_L_pillar_x-gap_horizontal-2.94,first_L_pillar_x-gap_horizontal-2.94],[first_L_pillar_y,first_L_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
           
       end
                
       %Line between corridors (Horizontal)
       plot([first_L_pillar_x,first_L_pillar_x-gap_horizontal],[first_L_pillar_y+(gap_vertical*12),first_L_pillar_y+(gap_vertical*12)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        
       plot([first_L_pillar_x,first_L_pillar_x-gap_horizontal],[first_L_pillar_y,first_L_pillar_y],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10); 
        
        %Upper corridor
        plot([first_L_pillar_x,first_L_pillar_x],[first_L_pillar_y+(gap_vertical*13),first_L_pillar_y+(gap_vertical*13)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        plot([first_R_pillar_x,first_R_pillar_x],[first_L_pillar_y+(gap_vertical*12),first_L_pillar_y+(gap_vertical*13)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10); 
        plot([first_L_pillar_x,first_L_pillar_x-gap_horizontal],[first_L_pillar_y+(gap_vertical*13),first_L_pillar_y+(gap_vertical*13)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        plot([first_L_pillar_x-gap_horizontal-2.94,first_L_pillar_x-gap_horizontal-2.94],[first_L_pillar_y+(gap_vertical*12),first_L_pillar_y+(gap_vertical*13)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        
        
    end
    
    if strcmp(type,'start_tag_28')    %Straight cw start from big tag 
        
%         plot(first_R_pillar_x,first_R_pillar_y,'--rs','LineWidth',2,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','g',...
%                     'MarkerSize',10);
%        plot(first_L_pillar_x,first_L_pillar_y,'--rs','LineWidth',2,...
%             'MarkerEdgeColor','k',...
%             'MarkerFaceColor','g',...
%             'MarkerSize',10);

       %Corridor A
       
       %axis([-2 7 -5 42]);
       
       for i=1:12
           
           plot([first_R_pillar_x,first_R_pillar_x],[first_R_pillar_y,first_R_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
           plot([first_L_pillar_x,first_L_pillar_x],[first_L_pillar_y,first_L_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
           
       end
       
       %Corridor B
       for i=1:12
           
           plot([first_R_pillar_x+gap_horizontal,first_R_pillar_x+gap_horizontal],[first_R_pillar_y,first_R_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
           plot([first_R_pillar_x+gap_horizontal+2.94,first_R_pillar_x+gap_horizontal+2.94],[first_L_pillar_y,first_L_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
           
       end
                
       %Line between corridors (Horizontal)
       plot([first_R_pillar_x,first_R_pillar_x+gap_horizontal],[first_L_pillar_y+(gap_vertical*12),first_L_pillar_y+(gap_vertical*12)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        
       plot([first_R_pillar_x,first_R_pillar_x+gap_horizontal],[first_L_pillar_y,first_L_pillar_y],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10); 
        
        %Upper corridor
        plot([first_L_pillar_x,first_L_pillar_x],[first_L_pillar_y+(gap_vertical*12),first_L_pillar_y+(gap_vertical*13)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        plot([first_R_pillar_x,first_R_pillar_x],[first_L_pillar_y+(gap_vertical*13),first_L_pillar_y+(gap_vertical*13)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10); 
        plot([first_R_pillar_x,first_R_pillar_x+gap_horizontal],[first_L_pillar_y+(gap_vertical*13),first_L_pillar_y+(gap_vertical*13)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        plot([first_R_pillar_x+gap_horizontal+2.94,first_R_pillar_x+gap_horizontal+2.94],[first_L_pillar_y+(gap_vertical*12),first_L_pillar_y+(gap_vertical*13)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        
        
    end
    
    if strcmp(type,'start_tag_22')    %Straight cw start from big tag 
        
%         axis([-7 2 -5 42]);
       
       for i=1:7
           
           plot([first_R_pillar_x,first_R_pillar_x],[first_R_pillar_y,first_R_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
           plot([first_L_pillar_x,first_L_pillar_x],[first_L_pillar_y,first_L_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
           
       end
       
       %Corridor A
       for i=1:7
           
           plot([first_L_pillar_x-gap_horizontal,first_L_pillar_x-gap_horizontal],[first_R_pillar_y,first_R_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
           plot([first_L_pillar_x-gap_horizontal-2.94,first_L_pillar_x-gap_horizontal-2.94],[first_L_pillar_y,first_L_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
           
       end
                
       %Line between corridors (Horizontal)
       plot([first_L_pillar_x,first_L_pillar_x-gap_horizontal],[first_L_pillar_y+(gap_vertical*7),first_L_pillar_y+(gap_vertical*7)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        
    end
    
    
    if strcmp(type,'start_tag_21')    %Straight cw start from big tag 
        
%         axis([-7 2 -5 42]);
       
       for i=1:8
           
           plot([first_R_pillar_x,first_R_pillar_x],[first_R_pillar_y,first_R_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
           plot([first_L_pillar_x,first_L_pillar_x],[first_L_pillar_y,first_L_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
           
       end
       
       %Corridor A
       for i=1:8
           
           plot([first_L_pillar_x-gap_horizontal,first_L_pillar_x-gap_horizontal],[first_R_pillar_y,first_R_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
           plot([first_L_pillar_x-gap_horizontal-2.94,first_L_pillar_x-gap_horizontal-2.94],[first_L_pillar_y,first_L_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
           
       end
                
       %Line between corridors (Horizontal)
       plot([first_L_pillar_x,first_L_pillar_x-gap_horizontal],[first_L_pillar_y+(gap_vertical*8),first_L_pillar_y+(gap_vertical*8)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        
    end
    
    if strcmp(type,'start_tag_16')    %Straight ccw start from big tag 
        
%         plot(first_R_pillar_x,first_R_pillar_y,'--rs','LineWidth',2,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','g',...
%                     'MarkerSize',10);
%        plot(first_L_pillar_x,first_L_pillar_y,'--rs','LineWidth',2,...
%             'MarkerEdgeColor','k',...
%             'MarkerFaceColor','g',...
%             'MarkerSize',10);

       %Corridor B 
       
%        axis([-7 2 -5 42]);
       
       for i=1:12
           
           plot([first_R_pillar_x,first_R_pillar_x],[first_R_pillar_y,first_R_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
           plot([first_L_pillar_x,first_L_pillar_x],[first_L_pillar_y,first_L_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
           
       end
       
       %Corridor A
       for i=1:12
           
           plot([first_L_pillar_x-gap_horizontal,first_L_pillar_x-gap_horizontal],[first_R_pillar_y,first_R_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',10);
                
           plot([first_L_pillar_x-gap_horizontal-2.94,first_L_pillar_x-gap_horizontal-2.94],[first_L_pillar_y,first_L_pillar_y+(gap_vertical*i)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
           
       end
                
       %Line between corridors (Horizontal)
       plot([first_L_pillar_x,first_L_pillar_x-gap_horizontal],[first_L_pillar_y+(gap_vertical*12),first_L_pillar_y+(gap_vertical*12)],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10);
        
       plot([first_L_pillar_x,first_L_pillar_x-gap_horizontal],[first_L_pillar_y,first_L_pillar_y],'--rs','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',10); 
        
        
        
    end
    
end