function plotCurrentData(app)
            subplot(1,3,1);
%            plot(time,Xdata,'b')
            timeline = app.fullTimeline_;
           plot(timeline, app.Zsim_, '.');
            title ('Z Position and Full Extension');
            xlabel ('Time (sec)');
            ylabel ('Bead Position (nm)');
end