
for i = 1:4
    output = Evaluate(i);
    
    figure('Color', 'w')
    subplot(2, 2, 1)
    plot(output.x_cur)
    title('x cur')
    
    subplot(2, 2, 2)
    plot(output.x_dstbd)
    title('x dstbd')
    
    subplot(2, 2, 3)
    plot(output.y); hold on;
    plot(output.y_d);
    title('y y_d')
    
    subplot(2, 2, 4)
    plot(output.I); hold on;
    title('I')
    
end