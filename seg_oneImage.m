function seg_oneImage(filename, minR, maxR)

	% segmentation of only one image

	f_path_trans     = filename;

	% initialize HoughTracker class
	% argument for the constructor is a name-string
	% output images are saved in the folder "name-string"_tracker
	ht = HoughTracker('example','cells_radii', [minR, maxR]);

	% rad the trans image
	img_trans = uint16( imread( f_path_trans ) );
	% do another time point
	% ht.add_timepoint( img_trans, 0, 0, 0 );
	[accum, circen, cirrad] = ht.hough_transform( img_trans );

	% write out circen and cirrad to file
	fileID = fopen( '/tmp/example_rois.txt','w' );
	formatSpec = '%4.4f\t%4.4f\t%3.1f\n'; % format?
	for i=1:size(circen,1)
		fprintf( fileID, formatSpec, circen(i,1), circen(i,2), cirrad(i) );
	end
	fclose( fileID );
end
