function rr = drawbullseye(cur_endo, cur_epi, difference, caseName)
    rr = 0;
    figlist = [131,132,133];
    close;
    figure('pos',[10 10 1500 350])
    for t = 133%figlist
        t
        switch t
            case 131
                coo = 'black';
                cur = cur_endo; % rescale(cur_endo,-1,1);
                cur_val = cur .* (abs(cur) >= 1.96.*stdd(1:17));
            case 132
                coo = 'black';
                cur = cur_epi; % rescale(cur_endo,-1,1);
                cur_val = cur .* (abs(cur) >= 1.96.*stdd(18:end));
            case 133
                coo = 'black';
                cur = difference; % rescale(cur_endo,-1,1);
                cur_val = cur;
        end
            subplot(t)
            c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
            set(c,'Color','k','LineWidth',1)

            area16 = cur_val(1:6);
            fillBullseye(area16, 1.5,2,60,420);

            area712 = cur_val(7:12);
            fillBullseye(area712, 1,1.5,60,420);

            area1316 = cur_val(13:16);
            fillBullseye(area1316,0.5,1,45,405);

            area17 = cur_val(17);
            fillBullseye(area17,0,0.5,0,360);
            text(-0.2,0,num2str(cur(17),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(-0.3,-0.7,num2str(cur(15),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(-0.2,0.7,num2str(cur(13),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(-0.2,1.2,num2str(cur(7),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(-0.25,1.7,num2str(cur(1),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(-0.25,-1.2,num2str(cur(10),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(-0.2,-1.7,num2str(cur(4),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');

            text(-1,0,num2str(cur(14),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(0.53,0,num2str(cur(16),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');

            text(-1.3,0.7,num2str(cur(8),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(-1.3,-0.7,num2str(cur(9),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(0.9,0.5,num2str(cur(12),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(0.9,-0.5,num2str(cur(11),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');

            text(1.2,1,num2str(cur(6),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(1.4,-1,num2str(cur(5),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(-1.6,1,num2str(cur(2),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            text(-1.6,-1,num2str(cur(3),2),'Color',coo,'FontSize',8,'FontName','Calibri','FontWeight','Bold');
            uistack(c,'top');
            
            switch t
                case 131
                    title([caseName,' Endo --- Blue:+  Yellow:-']);
                case 132
                    title([caseName,' Epi --- Blue:+  Yellow:-']);
                case 133
                    title(' Expand/Shrink --- Blue: thick  Yellow: thin');
            end  
    end
end