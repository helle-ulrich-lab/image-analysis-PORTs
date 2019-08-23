# This script tracks foci in 2D over time.
# authors: Ronald Wong and Nicola Zilio

import sys
from os import listdir
from os.path import isfile, join
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from ij import IJ, ImagePlus, WindowManager
import fiji.plugin.trackmate.Settings as Settings
import fiji.plugin.trackmate.Model as Model
import fiji.plugin.trackmate.SelectionModel as SelectionModel
import fiji.plugin.trackmate.TrackMate as TrackMate
import fiji.plugin.trackmate.Logger as Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
import fiji.plugin.trackmate.tracking.sparselap.SparseLAPTrackerFactory as SparseLAPTrackerFactory
import fiji.plugin.trackmate.tracking.LAPUtils as LAPUtils
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
import fiji.plugin.trackmate.features.spot.SpotContrastAndSNRAnalyzerFactory as SpotContrastAndSNRAnalyzerFactory
import fiji.plugin.trackmate.features.spot.SpotContrastAndSNRAnalyzer as SpotContrastAndSNRAnalyzer
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzerFactory as SpotIntensityAnalyzerFactory
import fiji.plugin.trackmate.features.track.TrackSpeedStatisticsAnalyzer as TrackSpeedStatisticsAnalyzer
import fiji.plugin.trackmate.features.track.TrackBranchingAnalyzer as TrackBranchingAnalyzer
import fiji.plugin.trackmate.features.track.TrackDurationAnalyzer as TrackDurationAnalyzer
import fiji.plugin.trackmate.features.track.TrackSpotQualityFeatureAnalyzer as TrackSpotQualityFeatureAnalyzer

# Define path to list of movies
path_to_movies = "C:\\Users\\ronawong\\Desktop\\TEST\\data\\"
list_of_movies = [f for f in listdir(path_to_movies) if isfile(join(path_to_movies, f))]   

for movie in list_of_movies:
    
    id = path_to_movies + movie
    options = ImporterOptions()
    options.setId(id)
    options.setAutoscale(False)
    bf_imp = BF.openImagePlus(options)
    
    imp = bf_imp[0]  
    
    #-------------------------
    # Instantiate model object
    #-------------------------
       
    model = Model()
       
    # Set logger
    model.setLogger(Logger.IJ_LOGGER)
       
    #------------------------
    # Prepare settings object
    #------------------------
          
    settings = Settings()
    settings.setFrom(imp)
          
    # Configure detector
    settings.detectorFactory = LogDetectorFactory()
    settings.detectorSettings = {
        'DO_SUBPIXEL_LOCALIZATION' : True,
        'RADIUS' : 0.25,
        'TARGET_CHANNEL' : 1,
        'THRESHOLD' : 50.,
        'DO_MEDIAN_FILTERING' : True,
    } 
       
    # Configure tracker
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap()
    settings.trackerSettings['LINKING_MAX_DISTANCE'] = 0.6
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 0.6
    settings.trackerSettings['MAX_FRAME_GAP'] = 1
    settings.trackerSettings['ALLOW_TRACK_MERGING'] = True
    settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = True
    settings.trackerSettings['MERGING_MAX_DISTANCE'] = 0.2
    settings.trackerSettings['SPLITTING_MAX_DISTANCE'] = 0.2

	# Use a value of 1-2 work the best in most cases
    settings.trackerSettings['LINKING_FEATURE_PENALTIES'] = {"QUALITY":1.0}
	
    # Add the analyzers for some spot features
    settings.addSpotAnalyzerFactory(SpotIntensityAnalyzerFactory())
    settings.addSpotAnalyzerFactory(SpotContrastAndSNRAnalyzerFactory())
       
    # Add an analyzer for some track features, such as the track mean speed.
    settings.addTrackAnalyzer(TrackSpeedStatisticsAnalyzer())

    settings.addTrackAnalyzer(TrackBranchingAnalyzer())
    settings.addTrackAnalyzer(TrackDurationAnalyzer())
    settings.addTrackAnalyzer(TrackSpotQualityFeatureAnalyzer())
         
    # Use a value between 30-35 works the best
    track_filter = FeatureFilter('TRACK_MEAN_QUALITY', 100, True)
    settings.addTrackFilter(track_filter)
          
    settings.initialSpotFilterValue = 1
       
    # print(str(settings))
          
    #----------------------
    # Instantiate trackmate
    #----------------------
       
    trackmate = TrackMate(model, settings)
          
    #------------
    # Execute all
    #------------
       
    ok = trackmate.checkInput()
    if not ok:
	
        sys.exit(str(trackmate.getErrorMessage()))
        
    ok = trackmate.process()
    if not ok:
        fhandle = open(path_to_movies + movie[:-4] + '_spot_statistics.txt', 'w')
        fhandle.write("TRACK_ID\tSPOT_ID\tPOSITION_X\tPOSITION_Y\tPOSITION_Z\tFRAME\tQUALITY\tSNR\tMEAN_INSTENSITY\n")

        ftrackhandle = open(path_to_movies + movie[:-4] + '_track_statistics.txt', 'w')
        ftrackhandle.write("Label\tNUMBER_SPLITS\tNUMBER_MERGES\tTRACK_DURATION\n")
        continue
        #sys.exit(str(trackmate.getErrorMessage()))
          
    #----------------
    # Display results
    #----------------
       
    model.getLogger().log('Found ' + str(model.getTrackModel().nTracks(True)) + ' tracks.')
        
    selectionModel = SelectionModel(model)
    displayer =  HyperStackDisplayer(model, selectionModel, imp)
    displayer.render()
    displayer.refresh()

    # Create output files
    fhandle = open(path_to_movies + movie[:-4] + '_spot_statistics.txt', 'w')
    fhandle.write("TRACK_ID\tSPOT_ID\tPOSITION_X\tPOSITION_Y\tPOSITION_Z\tFRAME\tQUALITY\tSNR\tMEAN_INSTENSITY\n")

    ftrackhandle = open(path_to_movies + movie[:-4] + '_track_statistics.txt', 'w')
    ftrackhandle.write("Label\tNUMBER_SPLITS\tNUMBER_MERGES\tTRACK_DURATION\n")
       
    # The feature model, that stores edge and track features
    fm = model.getFeatureModel()
       
    for id in model.getTrackModel().trackIDs(True):
       
        # Fetch the track feature from the feature model
        v = fm.getTrackFeature(id, 'TRACK_MEAN_SPEED')
        d = fm.getTrackFeature(id, 'TRACK_DURATION')
        num_splits = fm.getTrackFeature(id, TrackBranchingAnalyzer.NUMBER_SPLITS)
        num_merges = fm.getTrackFeature(id, TrackBranchingAnalyzer.NUMBER_MERGES)
        model.getLogger().log('')
        model.getLogger().log("Image: " + path_to_movies + movie + ", Track: " + str(id))
        ftrackhandle.write(str(id) + "\t" + str(int(num_splits)) + "\t" + str(int(num_merges)) + "\t" + str(d) + "\n")
           
        track = model.getTrackModel().trackSpots(id)
        
        for spot in track:
            sid = spot.ID()
            
            # Fetch spot features directly from spot
            x=spot.getFeature('POSITION_X')
            y=spot.getFeature('POSITION_Y')
            z=spot.getFeature('POSITION_Z')
            t=spot.getFeature('FRAME')
            q=spot.getFeature('QUALITY')
            snr=spot.getFeature('SNR') 
            mean=spot.getFeature('MEAN_INTENSITY')
            fhandle.write(str(id) + "\t" + str(sid) + "\t" + str(x) + "\t" +str(y) + "\t" +str(z) + "\t" + str(int(t)) + "\t" + str(q) + "\t" + str(snr) + "\t" + str(mean) + "\n")

    fhandle.close()
    ftrackhandle.close()
    imp.close()