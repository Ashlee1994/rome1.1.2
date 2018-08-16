/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef METADATA_H_
#define METADATA_H_

#include "./util.h"
#include "./memory.h"
#include "./error.h"
#include "./string.h"
#include "./macros.h"
#include "./mpi.h"

// metadata name format number@filename ,
// like : 3@particles.mrcs
struct MetaName{
public:
    size_t INDEX;
    std::string NAME;
    MetaName():INDEX(0),NAME(""){}
    MetaName(const char* full_name){set(full_name);}
    MetaName(const std::string& full_name){set(full_name);}
    ~MetaName(){}
    bool operator!=(const MetaName& rhs) const {return !(rhs.INDEX==INDEX && rhs.NAME==NAME);}
    bool operator==(const MetaName& rhs) const {return  (rhs.INDEX==INDEX && rhs.NAME==NAME);}
    bool operator>(const MetaName& rhs) const {return  NAME > rhs.NAME ? true : INDEX > rhs.INDEX;}
    bool operator<(const MetaName& rhs) const {return  NAME < rhs.NAME ? true : INDEX < rhs.INDEX;}
private:
    void set(const std::string& full_name){
        INDEX = atoi(full_name.substr(0,full_name.find_last_of('@')).c_str());
        NAME = full_name.substr(full_name.find_last_of('@')+1);
    }
};
inline std::ostream& operator<<(std::ostream& os, MetaName const & rhs) {os<<num2str(rhs.INDEX,7)<<'@'<<rhs.NAME; return os; }

//  metadata elements
//          datatype   , var                , default   ,  enum (star file's header)
#define METADATA_ELTS \
    ELTDOU(FDOUBLE     , CTF_VOLTAGE		, 	0.0     ,   Voltage                     ) SEP \
    ELTDOU(FDOUBLE     , CTF_DEFOCUS_U		,   0.0     ,   DefocusU                    ) SEP \
    ELTDOU(FDOUBLE     , CTF_DEFOCUS_V		,   0.0     ,   DefocusV                    ) SEP \
    ELTDOU(FDOUBLE     , CTF_DEFOCUS_ANGLE	,   0.0     ,   DefocusAngle                ) SEP \
    ELTDOU(FDOUBLE     , CTF_CS             ,   0.0     ,   SphericalAberration         ) SEP \
    ELTDOU(FDOUBLE     , CTF_Q0             ,   0.0     ,   AmplitudeContrast           ) SEP \
    ELTDOU(FDOUBLE     , CTF_BFAC			,   0.0     ,   Bfactor                     ) SEP \
    \
    ELTSTR(MetaName    , IMAGE              , ""        ,   ImageName                   ) SEP \
    ELTSTR(std::string , MICROGRAPH         , ""        ,   MicrographName              ) SEP \
    /* for each group a different noise spectrum and signal scale factor is estimated */      \
    /* similar defocus values often have similar scale intensity factors and noise */         \
    /* power spectra,so if some groups image number is less,you can group some Micrograph */  \
    /* only the group_name for reading */                                                     \
    ELTINT(int         , GROUP_NO			, 1         ,   GroupNumber                 ) SEP \
    ELTSTR(std::string , GROUP_NAME			, ""        ,   GroupName                   ) SEP \
	/* Random subset this particle belongs to(only used in 'autorefine') */					  \
	ELTINT(int		   , SUBSET				, 0			,	RandomSubset				) SEP \
    ELTDOU(FDOUBLE     , ROT				, 0.0       ,   AngleRot                    ) SEP \
    ELTDOU(FDOUBLE     , TILT				, 0.0       ,   AngleTilt                   ) SEP \
    ELTDOU(FDOUBLE     , PSI				, 0.0       ,   AnglePsi                    ) SEP \
    ELTDOU(FDOUBLE     , XOFF				, 0.0       ,   OriginX                     ) SEP \
    ELTDOU(FDOUBLE     , YOFF				, 0.0       ,   OriginY                     ) SEP \
    ELTDOU(FDOUBLE     , ZOFF				, 0.0       ,   OriginZ                     ) SEP \
    ELTINT(int         , CLASS				, 0         ,   ClassNumber                 ) SEP \
    ELTDOU(FDOUBLE     , NORM				, 1.0       ,   NormCorrection              ) SEP \
    ELTDOU(FDOUBLE     , SCALE				, 1.0       ,   ScaleCorrection             ) SEP \
    ELTDOU(FDOUBLE     , DLL				, 0.0       ,   LogLikeliContribution       ) SEP \
    ELTDOU(FDOUBLE     , PMAX				, 0.0       ,   MaxValueProbDistribution    ) SEP \
    ELTINT(int         , NR_SIGN			, 0         ,   NrOfSignificantSamples      ) SEP \
    \
    ELTDOU(FDOUBLE     , ROT_PRIOR			, 999.0     ,   AngleRotPrior               ) SEP \
    ELTDOU(FDOUBLE     , TILT_PRIOR         , 999.0     ,   AngleTiltPrior              ) SEP \
    ELTDOU(FDOUBLE     , PSI_PRIOR			, 999.0     ,   AnglePsiPrior               ) SEP \
    ELTDOU(FDOUBLE     , XOFF_PRIOR         , 999.0     ,   OriginXPrior                ) SEP \
    ELTDOU(FDOUBLE     , YOFF_PRIOR         , 999.0     ,   OriginYPrior                ) SEP \
    ELTDOU(FDOUBLE     , ZOFF_PRIOR         , 999.0     ,   OriginZPrior                ) SEP \
    \
    ELTDOU(FDOUBLE     , MAGNIFICATION      ,   0.0     ,   Magnification               ) SEP \
    ELTDOU(FDOUBLE     , DETECTOR_PIXEL_SIZE,   0.0     ,   DetectorPixelSize           ) SEP \
    ELTDOU(FDOUBLE     , BEAMTILT_X         ,   0.0     ,   BeamTiltX                   ) SEP \
    ELTDOU(FDOUBLE     , BEAMTILT_Y         ,   0.0     ,   BeamTiltY                   ) SEP \
    \
    ELTDOU(FDOUBLE     , MAT_0_0			,   1.0     ,   Matrix_1_1                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_0_1			,   0.0     ,   Matrix_1_2                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_0_2			,   0.0     ,   Matrix_1_3                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_1_0			,   0.0     ,   Matrix_2_1                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_1_1			,   1.0     ,   Matrix_2_2                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_1_2			,   0.0     ,   Matrix_2_3                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_2_0			,   0.0     ,   Matrix_3_1                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_2_1			,   0.0     ,   Matrix_3_2                  ) SEP \
    ELTDOU(FDOUBLE     , MAT_2_2			,   1.0     ,   Matrix_3_3                  ) SEP \
    ELTSTR(std::string , AREANAME           ,   ""      ,   AreaName                    ) SEP \
    ELTSTR(std::string , BODYMASKNAME       ,   ""      , BodyMaskName                 ) SEP \
    ELTINT(int         , BODYKEEPFIXED      , 0         , BodyKeepFixed                ) SEP \
    ELTSTR(std::string , BODYREFERENCENAME  , ""        , BodyReferenceName            ) SEP \
    ELTDOU(FDOUBLE     , BODYROTATEDIRECTIONX, 0.0       , BodyRotateDirectionX         ) SEP \
    ELTDOU(FDOUBLE     , BODYROTATEDIRECTIONY, 0.0       , BodyRotateDirectionY         ) SEP \
    ELTDOU(FDOUBLE     , BODYROTATEDIRECTIONZ, 0.0       , BodyRotateDirectionZ         ) SEP \
    ELTINT(int         , BODYROTATERELATIVETO, 0         , BodyRotateRelativeTo         ) SEP \
    ELTDOU(FDOUBLE     , BODYSIGMAANGLES    , 0.0       , BodySigmaAngles              ) SEP \
    ELTDOU(FDOUBLE     , BODYSIGMAOFFSET    , 0.0       , BodySigmaOffset              ) SEP \
    ELTDOU(FDOUBLE     , BODYSIGMAROT       , 0.0       , BodySigmaRot                 ) SEP \
    ELTDOU(FDOUBLE     , BODYSIGMATILT      , 0.0       , BodySigmaTilt                ) SEP \
    ELTDOU(FDOUBLE     , BODYSIGMAPSI       , 0.0       , BodySigmaPsi                 ) SEP \
    ELTSTR(std::string , BODYSTARFILE       , ""        , BodyStarFile                 ) SEP \
    ELTDOU(FDOUBLE     , CTFASTIGMATISM     , 0.0       , CtfAstigmatism               ) SEP \
    ELTDOU(FDOUBLE     , CTFBFACTOR         , 0.0       , CtfBfactor                   ) SEP \
    ELTDOU(FDOUBLE     , CTFMAXRESOLUTION   , 0.0       , CtfMaxResolution             ) SEP \
    ELTDOU(FDOUBLE     , CTFVALIDATIONSCORE , 0.0       , CtfValidationScore           ) SEP \
    ELTDOU(FDOUBLE     , CTFSCALEFACTOR     , 0.0       , CtfScalefactor               ) SEP \
    ELTDOU(FDOUBLE     , CHROMATICABERRATION, 0.0       , ChromaticAberration          ) SEP \
    ELTDOU(FDOUBLE     , ENERGYLOSS         , 0.0       , EnergyLoss                   ) SEP \
    ELTDOU(FDOUBLE     , CTFFIGUREOFMERIT   , 0.0       , CtfFigureOfMerit             ) SEP \
    ELTSTR(std::string , CTFIMAGE           , ""        , CtfImage                     ) SEP \
    ELTDOU(FDOUBLE     , LENSSTABILITY      , 0.0       , LensStability                ) SEP \
    ELTDOU(FDOUBLE     , PHASESHIFT         , 0.0       , PhaseShift                   ) SEP \
    ELTDOU(FDOUBLE     , CONVERGENCECONE    , 0.0       , ConvergenceCone              ) SEP \
    ELTDOU(FDOUBLE     , LONGITUDINALDISPLACEMENT, 0.0  , LongitudinalDisplacement     ) SEP \
    ELTDOU(FDOUBLE     , TRANSVERSALDISPLACEMENT, 0.0   , TransversalDisplacement      ) SEP \
    ELTDOU(FDOUBLE     , CTFVALUE           , 0.0       , CtfValue                     ) SEP \
    ELTSTR(std::string , IMAGEORIGINALNAME  , ""        , ImageOriginalName            ) SEP \
    ELTSTR(std::string , RECONSTRUCTIMAGENAME, ""       , ReconstructImageName         ) SEP \
    ELTINT(int         , IMAGEID            , 0         , ImageId                      ) SEP \
    ELTINT(int         , ENABLED            , 1         , Enabled                      ) SEP \
    ELTINT(int         , DATATYPE           , 0         , DataType                     ) SEP \
    ELTINT(int         , IMAGEDIMENSIONALITY, 0         , ImageDimensionality          ) SEP \
    ELTSTR(std::string , BEAMTILTGROUPNAME  , ""        , BeamTiltGroupName            ) SEP \
    ELTDOU(FDOUBLE     , COORDINATEX        , 0.0       , CoordinateX                  ) SEP \
    ELTDOU(FDOUBLE     , COORDINATEY        , 0.0       , CoordinateY                  ) SEP \
    ELTDOU(FDOUBLE     , COORDINATEZ        , 0.0       , CoordinateZ                  ) SEP \
    ELTINT(int         , MOVIEFRAMENUMBER   , 0         , MovieFrameNumber             ) SEP \
    ELTDOU(FDOUBLE     , MAGNIFICATIONCORRECTION, 0.0       , MagnificationCorrection      ) SEP \
    ELTDOU(FDOUBLE     , SAMPLINGRATE       , 0.0       , SamplingRate                 ) SEP \
    ELTDOU(FDOUBLE     , SAMPLINGRATEX      , 0.0       , SamplingRateX                ) SEP \
    ELTDOU(FDOUBLE     , SAMPLINGRATEY      , 0.0       , SamplingRateY                ) SEP \
    ELTDOU(FDOUBLE     , SAMPLINGRATEZ      , 0.0       , SamplingRateZ                ) SEP \
    ELTINT(int         , IMAGESIZE          , 0         , ImageSize                    ) SEP \
    ELTINT(int         , IMAGESIZEX         , 0         , ImageSizeX                   ) SEP \
    ELTINT(int         , IMAGESIZEY         , 0         , ImageSizeY                   ) SEP \
    ELTINT(int         , IMAGESIZEZ         , 0         , ImageSizeZ                   ) SEP \
    ELTDOU(FDOUBLE     , MINIMUMVALUE       , 0.0       , MinimumValue                 ) SEP \
    ELTDOU(FDOUBLE     , MAXIMUMVALUE       , 0.0       , MaximumValue                 ) SEP \
    ELTDOU(FDOUBLE     , AVERAGEVALUE       , 0.0       , AverageValue                 ) SEP \
    ELTDOU(FDOUBLE     , STANDARDDEVIATIONVALUE, 0.0       , StandardDeviationValue       ) SEP \
    ELTDOU(FDOUBLE     , SKEWNESSVALUE      , 0.0       , SkewnessValue                ) SEP \
    ELTDOU(FDOUBLE     , KURTOSISEXCESSVALUE, 0.0       , KurtosisExcessValue          ) SEP \
    ELTDOU(FDOUBLE     , IMAGEWEIGHT        , 0.0       , ImageWeight                  ) SEP \
    ELTSTR(std::string , MASKNAME           , ""        , MaskName                     ) SEP \
    ELTDOU(FDOUBLE     , ACCUMMOTIONTOTAL   , 0.0       , AccumMotionTotal             ) SEP \
    ELTDOU(FDOUBLE     , ACCUMMOTIONEARLY   , 0.0       , AccumMotionEarly             ) SEP \
    ELTDOU(FDOUBLE     , ACCUMMOTIONLATE    , 0.0       , AccumMotionLate              ) SEP \
    ELTINT(int         , MICROGRAPHID       , 0         , MicrographId                 ) SEP \
    ELTSTR(std::string , MICROGRAPHGAINNAME , ""        , MicrographGainName           ) SEP \
    ELTSTR(std::string , MICROGRAPHDEFECTFILE, ""        , MicrographDefectFile         ) SEP \
    ELTSTR(std::string , MICROGRAPHNAMENODW , ""        , MicrographNameNoDW           ) SEP \
    ELTSTR(std::string , MICROGRAPHMOVIENAME, ""        , MicrographMovieName          ) SEP \
    ELTSTR(std::string , MICROGRAPHMETADATA , ""        , MicrographMetadata           ) SEP \
    ELTDOU(FDOUBLE     , MICROGRAPHTILTANGLE, 0.0       , MicrographTiltAngle          ) SEP \
    ELTDOU(FDOUBLE     , MICROGRAPHTILTAXISDIRECTION, 0.0       , MicrographTiltAxisDirection  ) SEP \
    ELTDOU(FDOUBLE     , MICROGRAPHTILTAXISOUTOFPLANE, 0.0       , MicrographTiltAxisOutOfPlane ) SEP \
    ELTDOU(FDOUBLE     , MICROGRAPHORIGINALPIXELSIZE, 0.0       , MicrographOriginalPixelSize  ) SEP \
    ELTDOU(FDOUBLE     , MICROGRAPHPREEXPOSURE, 0.0       , MicrographPreExposure        ) SEP \
    ELTDOU(FDOUBLE     , MICROGRAPHDOSERATE , 0.0       , MicrographDoseRate           ) SEP \
    ELTDOU(FDOUBLE     , MICROGRAPHBINNING  , 0.0       , MicrographBinning            ) SEP \
    ELTINT(int         , MICROGRAPHFRAMENUMBER, 0         , MicrographFrameNumber        ) SEP \
    ELTINT(int         , MOTIONMODELVERSION , 0         , MotionModelVersion           ) SEP \
    ELTINT(int         , MICROGRAPHSTARTFRAME, 0         , MicrographStartFrame         ) SEP \
    ELTINT(int         , MICROGRAPHENDFRAME , 0         , MicrographEndFrame           ) SEP \
    ELTDOU(FDOUBLE     , MICROGRAPHSHIFTX   , 0.0       , MicrographShiftX             ) SEP \
    ELTDOU(FDOUBLE     , MICROGRAPHSHIFTY   , 0.0       , MicrographShiftY             ) SEP \
    ELTINT(int         , MOTIONMODELCOEFFSIDX, 0         , MotionModelCoeffsIdx         ) SEP \
    ELTDOU(FDOUBLE     , MOTIONMODELCOEFF   , 0.0       , MotionModelCoeff             ) SEP \
    ELTDOU(FDOUBLE     , ACCURACYROTATIONS  , 0.0       , AccuracyRotations            ) SEP \
    ELTDOU(FDOUBLE     , ACCURACYTRANSLATIONS, 0.0       , AccuracyTranslations         ) SEP \
    ELTDOU(FDOUBLE     , AVERAGEPMAX        , 0.0       , AveragePmax                  ) SEP \
    ELTDOU(FDOUBLE     , CURRENTRESOLUTION  , 0.0       , CurrentResolution            ) SEP \
    ELTINT(int         , CURRENTIMAGESIZE   , 0         , CurrentImageSize             ) SEP \
    ELTDOU(FDOUBLE     , SSNRMAP            , 0.0       , SsnrMap                      ) SEP \
    ELTINT(int         , REFERENCEDIMENSIONALITY, 0         , ReferenceDimensionality      ) SEP \
    ELTINT(int         , DATADIMENSIONALITY , 0         , DataDimensionality           ) SEP \
    ELTDOU(FDOUBLE     , DIFF2RANDOMHALVES  , 0.0       , Diff2RandomHalves            ) SEP \
    ELTDOU(FDOUBLE     , ESTIMATEDRESOLUTION, 0.0       , EstimatedResolution          ) SEP \
    ELTDOU(FDOUBLE     , FOURIERCOMPLETENESS, 0.0       , FourierCompleteness          ) SEP \
    ELTDOU(FDOUBLE     , OVERALLFOURIERCOMPLETENESS, 0.0       , OverallFourierCompleteness   ) SEP \
    ELTDOU(FDOUBLE     , GOLDSTANDARDFSC    , 0.0       , GoldStandardFsc              ) SEP \
    ELTINT(int         , GROUPNRPARTICLES   , 0         , GroupNrParticles             ) SEP \
    ELTDOU(FDOUBLE     , GROUPSCALECORRECTION, 0.0       , GroupScaleCorrection         ) SEP \
    ELTINT(int         , NRHELICALASYMUNITS , 0         , NrHelicalAsymUnits           ) SEP \
    ELTDOU(FDOUBLE     , HELICALTWIST       , 0.0       , HelicalTwist                 ) SEP \
    ELTDOU(FDOUBLE     , HELICALTWISTMIN    , 0.0       , HelicalTwistMin              ) SEP \
    ELTDOU(FDOUBLE     , HELICALTWISTMAX    , 0.0       , HelicalTwistMax              ) SEP \
    ELTDOU(FDOUBLE     , HELICALTWISTINITIALSTEP, 0.0       , HelicalTwistInitialStep      ) SEP \
    ELTDOU(FDOUBLE     , HELICALRISE        , 0.0       , HelicalRise                  ) SEP \
    ELTDOU(FDOUBLE     , HELICALRISEMIN     , 0.0       , HelicalRiseMin               ) SEP \
    ELTDOU(FDOUBLE     , HELICALRISEMAX     , 0.0       , HelicalRiseMax               ) SEP \
    ELTDOU(FDOUBLE     , HELICALRISEINITIALSTEP, 0.0       , HelicalRiseInitialStep       ) SEP \
    ELTINT(FDOUBLE     , ISHELIX            , 1         , IsHelix                      ) SEP \
    ELTINT(int         , FOURIERSPACEINTERPOLATOR, 0         , FourierSpaceInterpolator     ) SEP \
    ELTDOU(FDOUBLE     , LOGLIKELIHOOD      , 0.0       , LogLikelihood                ) SEP \
    ELTINT(int         , MINRADIUSNNINTERPOLATION, 0         , MinRadiusNnInterpolation ) SEP \
    ELTDOU(FDOUBLE     , NORMCORRECTIONAVERAGE, 0.0       , NormCorrectionAverage        ) SEP \
    ELTINT(int         , NRCLASSES          , 0         , NrClasses                    ) SEP \
    ELTINT(int         , NRBODIES           , 0         , NrBodies                     ) SEP \
    ELTINT(int         , NRGROUPS           , 0         , NrGroups                     ) SEP \
    ELTDOU(FDOUBLE     , SPECTRALORIENTABILITYCONTRIBUTION, 0.0       , SpectralOrientabilityContribution ) SEP \
    ELTINT(int         , ORIGINALIMAGESIZE  , 0         , OriginalImageSize            ) SEP \
    ELTDOU(FDOUBLE     , PADDINGFACTOR      , 0.0       , PaddingFactor                ) SEP \
    ELTDOU(FDOUBLE     , CLASSDISTRIBUTION  , 0.0       , ClassDistribution            ) SEP \
    ELTDOU(FDOUBLE     , CLASSPRIOROFFSETX  , 0.0       , ClassPriorOffsetX            ) SEP \
    ELTDOU(FDOUBLE     , CLASSPRIOROFFSETY  , 0.0       , ClassPriorOffsetY            ) SEP \
    ELTDOU(FDOUBLE     , ORIENTATIONDISTRIBUTION, 0.0       , OrientationDistribution      ) SEP \
    ELTDOU(FDOUBLE     , PIXELSIZE          , 0.0       , PixelSize                    ) SEP \
    ELTDOU(FDOUBLE     , REFERENCESPECTRALPOWER, 0.0       , ReferenceSpectralPower       ) SEP \
    ELTINT(int         , ORIENTATIONALPRIORMODE, 0         , OrientationalPriorMode       ) SEP \
    ELTSTR(std::string , REFERENCEIMAGE     , ""        , ReferenceImage               ) SEP \
    ELTSTR(std::string , SGDGRADIENTIMAGE   , ""        , SGDGradientImage             ) SEP \
    ELTDOU(FDOUBLE     , SIGMAOFFSETS       , 0.0       , SigmaOffsets                 ) SEP \
    ELTDOU(FDOUBLE     , SIGMA2NOISE        , 0.0       , Sigma2Noise                  ) SEP \
    ELTDOU(FDOUBLE     , REFERENCESIGMA2    , 0.0       , ReferenceSigma2              ) SEP \
    ELTDOU(FDOUBLE     , SIGMAPRIORROTANGLE , 0.0       , SigmaPriorRotAngle           ) SEP \
    ELTDOU(FDOUBLE     , SIGMAPRIORTILTANGLE, 0.0       , SigmaPriorTiltAngle          ) SEP \
    ELTDOU(FDOUBLE     , SIGMAPRIORPSIANGLE , 0.0       , SigmaPriorPsiAngle           ) SEP \
    ELTDOU(FDOUBLE     , SIGNALTONOISERATIO , 0.0       , SignalToNoiseRatio           ) SEP \
    ELTDOU(FDOUBLE     , TAU2FUDGEFACTOR    , 0.0       , Tau2FudgeFactor              ) SEP \
    ELTDOU(FDOUBLE     , REFERENCETAU2      , 0.0       , ReferenceTau2                ) SEP \
    ELTDOU(FDOUBLE     , OVERALLACCURACYROTATIONS, 0.0       , OverallAccuracyRotations     ) SEP \
    ELTDOU(FDOUBLE     , OVERALLACCURACYTRANSLATIONS, 0.0       , OverallAccuracyTranslations  ) SEP \
    ELTDOU(FDOUBLE     , ADAPTIVEOVERSAMPLEFRACTION, 0.0       , AdaptiveOversampleFraction   ) SEP \
    ELTINT(int         , ADAPTIVEOVERSAMPLEORDER, 0         , AdaptiveOversampleOrder      ) SEP \
    ELTINT(int         , AUTOLOCALSEARCHESHEALPIXORDER, 0         , AutoLocalSearchesHealpixOrder ) SEP \
    ELTDOU(FDOUBLE     , AVAILABLEMEMORY    , 0.0       , AvailableMemory              ) SEP \
    ELTDOU(FDOUBLE     , BESTRESOLUTIONTHUSFAR, 0.0       , BestResolutionThusFar        ) SEP \
    ELTINT(int         , COARSEIMAGESIZE    , 0         , CoarseImageSize              ) SEP \
    ELTDOU(FDOUBLE     , CHANGESOPTIMALOFFSETS, 0.0       , ChangesOptimalOffsets        ) SEP \
    ELTDOU(FDOUBLE     , CHANGESOPTIMALORIENTATIONS, 0.0       , ChangesOptimalOrientations   ) SEP \
	ELTDOU(FDOUBLE     , CHANGESOPTIMALCLASSES, 0.0       , ChangesOptimalClasses        ) SEP \
    ELTINT(FDOUBLE     , CTFDATAAREPHASEFLIPPED, 1         , CtfDataArePhaseFlipped       ) SEP \
    ELTINT(FDOUBLE     , CTFDATAARECTFPREMULTIPLIED, 1         , CtfDataAreCtfPremultiplied   ) SEP \
    ELTSTR(std::string , EXPERIMENTALDATASTARFILE, ""        , ExperimentalDataStarFile     ) SEP \
    ELTINT(bool , DOCORRECTCTF       , 1         , DoCorrectCtf                 ) SEP \
    ELTINT(bool , DOCORRECTMAGNIFICATION, 1         , DoCorrectMagnification       ) SEP \
    ELTINT(bool , DOCORRECTNORM      , 1         , DoCorrectNorm                ) SEP \
    ELTINT(bool , DOCORRECTSCALE     , 1         , DoCorrectScale               ) SEP \
	ELTINT(bool , DOREALIGNMOVIES    , 1         , DoRealignMovies              ) SEP \
    ELTINT(bool		   , DOMAPESTIMATION    , 1         , DoMapEstimation              ) SEP \
    ELTINT(bool 	   , DOSTOCHASTICGRADIENTDESCENT, 1         , DoStochasticGradientDescent  ) SEP \
    ELTINT(bool		   , DOFASTSUBSETOPTIMISATION, 1         , DoFastSubsetOptimisation     ) SEP \
    ELTINT(int         , SGDINITIALITERATIONS, 0         , SgdInitialIterations         ) SEP \
    ELTINT(int         , SGDFINALITERATIONS , 0         , SgdFinalIterations           ) SEP \
    ELTINT(int         , SGDINBETWEENITERATIONS, 0         , SgdInBetweenIterations       ) SEP \
    ELTDOU(FDOUBLE     , SGDINITIALRESOLUTION, 0.0       , SgdInitialResolution         ) SEP \
    ELTDOU(FDOUBLE     , SGDFINALRESOLUTION , 0.0       , SgdFinalResolution           ) SEP \
    ELTINT(int         , SGDINITIALSUBSETSIZE, 0         , SgdInitialSubsetSize         ) SEP \
    ELTINT(int         , SGDFINALSUBSETSIZE , 0         , SgdFinalSubsetSize           ) SEP \
    ELTDOU(FDOUBLE     , SGDMUFACTOR        , 0.0       , SgdMuFactor                  ) SEP \
    ELTDOU(FDOUBLE     , SGDSIGMA2FUDGEINITIAL, 0.0       , SgdSigma2FudgeInitial        ) SEP \
    ELTINT(int         , SGDSIGMA2FUDGEHALFLIFE, 0         , SgdSigma2FudgeHalflife       ) SEP \
    ELTINT(int         , SGDSKIPANNEAL      , 1         , SgdSkipAnneal                ) SEP \
    ELTINT(int         , SGDSUBSETSIZE      , 0         , SgdSubsetSize                ) SEP \
    ELTINT(int         , SGDWRITEEVERYSUBSET, 0         , SgdWriteEverySubset          ) SEP \
    ELTINT(int         , SGDMAXSUBSETS      , 0         , SgdMaxSubsets                ) SEP \
    ELTDOU(FDOUBLE     , SGDSTEPSIZE        , 0.0       , SgdStepsize                  ) SEP \
    ELTINT(FDOUBLE     , DOAUTOREFINE       , 1         , DoAutoRefine                 ) SEP \
    ELTINT(bool        , DOONLYFLIPCTFPHASES, 1         , DoOnlyFlipCtfPhases          ) SEP \
    ELTINT(bool        , DOSOLVENTFLATTENING, 1         , DoSolventFlattening          ) SEP \
    ELTINT(bool        , DOSOLVENTFSCCORRECTION, 1         , DoSolventFscCorrection       ) SEP \
    ELTINT(bool        , DOSKIPALIGN        , 1         , DoSkipAlign                  ) SEP \
    ELTINT(bool        , DOSKIPROTATE       , 1         , DoSkipRotate                 ) SEP \
    ELTINT(bool        , DOSPLITRANDOMHALVES, 1         , DoSplitRandomHalves          ) SEP \
    ELTINT(bool        , DOZEROMASK         , 1         , DoZeroMask                   ) SEP \
    ELTINT(bool        , FIXSIGMANOISEESTIMATES, 1         , FixSigmaNoiseEstimates       ) SEP \
    ELTINT(bool        , FIXSIGMAOFFSETESTIMATES, 1         , FixSigmaOffsetEstimates      ) SEP \
    ELTINT(bool        , FIXTAUESTIMATES    , 1         , FixTauEstimates              ) SEP \
    ELTINT(bool        , HASCONVERGED       , 1         , HasConverged                 ) SEP \
    ELTINT(bool        , HASHIGHFSCATRESOLLIMIT, 1         , HasHighFscAtResolLimit       ) SEP \
    ELTINT(int         , HASLARGESIZEINCREASEITERATIONSAGO, 0         , HasLargeSizeIncreaseIterationsAgo ) SEP \
    ELTINT(bool        , DOHELICALREFINE    , 1         , DoHelicalRefine              ) SEP \
    ELTINT(bool        , IGNOREHELICALSYMMETRY, 1         , IgnoreHelicalSymmetry        ) SEP \
    ELTDOU(FDOUBLE     , HELICALTWISTINITIAL, 0.0       , HelicalTwistInitial          ) SEP \
    ELTDOU(FDOUBLE     , HELICALRISEINITIAL , 0.0       , HelicalRiseInitial           ) SEP \
    ELTDOU(FDOUBLE     , HELICALCENTRALPROPORTION, 0.0       , HelicalCentralProportion     ) SEP \
    ELTDOU(FDOUBLE     , HELICALMASKTUBEINNERDIAMETER, 0.0       , HelicalMaskTubeInnerDiameter ) SEP \
    ELTDOU(FDOUBLE     , HELICALMASKTUBEOUTERDIAMETER, 0.0       , HelicalMaskTubeOuterDiameter ) SEP \
    ELTINT(bool        , HELICALSYMMETRYLOCALREFINEMENT, 1         , HelicalSymmetryLocalRefinement ) SEP \
    ELTDOU(FDOUBLE     , HELICALSIGMADISTANCE, 0.0       , HelicalSigmaDistance         ) SEP \
    ELTINT(bool        , HELICALKEEPTILTPRIORFIXED, 1         , HelicalKeepTiltPriorFixed    ) SEP \
    ELTDOU(FDOUBLE     , HIGHRESLIMITEXPECTATION, 0.0       , HighresLimitExpectation      ) SEP \
    ELTDOU(FDOUBLE     , HIGHRESLIMITSGD    , 0.0       , HighresLimitSGD              ) SEP \
    ELTINT(bool        , DOIGNORECTFUNTILFIRSTPEAK, 1         , DoIgnoreCtfUntilFirstPeak    ) SEP \
    ELTINT(int         , INCREMENTIMAGESIZE , 0         , IncrementImageSize           ) SEP \
    ELTINT(int         , CURRENTITERATION   , 0         , CurrentIteration             ) SEP \
    ELTSTR(std::string , LOCALSYMMETRYFILE  , ""        , LocalSymmetryFile            ) SEP \
    ELTDOU(FDOUBLE     , JOINHALVESUNTILTHISRESOLUTION, 0.0       , JoinHalvesUntilThisResolution ) SEP \
    ELTDOU(FDOUBLE     , MAGNIFICATIONSEARCHRANGE, 0.0       , MagnificationSearchRange     ) SEP \
    ELTDOU(FDOUBLE     , MAGNIFICATIONSEARCHSTEP, 0.0       , MagnificationSearchStep      ) SEP \
    ELTINT(int         , MAXIMUMCOARSEIMAGESIZE, 0         , MaximumCoarseImageSize       ) SEP \
    ELTINT(int         , MAXNUMBEROFPOOLEDPARTICLES, 0         , MaxNumberOfPooledParticles   ) SEP \
    ELTSTR(std::string , MODELSTARFILE      , ""        , ModelStarFile                ) SEP \
    ELTSTR(std::string , MODELSTARFILE2     , ""        , ModelStarFile2               ) SEP \
    ELTINT(int         , NUMBEROFITERATIONS , 0         , NumberOfIterations           ) SEP \
    ELTINT(int         , NUMBEROFITERWITHOUTRESOLUTIONGAIN, 0         , NumberOfIterWithoutResolutionGain ) SEP \
    ELTINT(int         , NUMBEROFITERWITHOUTCHANGINGASSIGNMENTS, 0         , NumberOfIterWithoutChangingAssignments ) SEP \
    ELTSTR(std::string , OUTPUTROOTNAME     , ""        , OutputRootName               ) SEP \
    ELTDOU(FDOUBLE     , PARTICLEDIAMETER   , 0.0       , ParticleDiameter             ) SEP \
    ELTINT(int         , RADIUSMASKMAP      , 0         , RadiusMaskMap                ) SEP \
    ELTINT(int         , RADIUSMASKEXPIMAGES, 0         , RadiusMaskExpImages          ) SEP \
    ELTINT(int         , RANDOMSEED         , 0         , RandomSeed                   ) SEP \
    ELTINT(bool        , REFSARECTFCORRECTED, 1         , RefsAreCtfCorrected          ) SEP \
    ELTINT(int         , SMALLESTCHANGESCLASSES, 0         , SmallestChangesClasses       ) SEP \
    ELTDOU(FDOUBLE     , SMALLESTCHANGESOFFSETS, 0.0       , SmallestChangesOffsets       ) SEP \
    ELTDOU(FDOUBLE     , SMALLESTCHANGESORIENTATIONS, 0.0       , SmallestChangesOrientations  ) SEP \
    ELTSTR(std::string , ORIENTSAMPLINGSTARFILE, ""        , OrientSamplingStarFile       ) SEP \
    ELTSTR(std::string , SOLVENTMASKNAME    , ""        , SolventMaskName              ) SEP \
    ELTSTR(std::string , SOLVENTMASK2NAME   , ""        , SolventMask2Name             ) SEP \
    ELTSTR(std::string , TAUSPECTRUMNAME    , ""        , TauSpectrumName              ) SEP \
    ELTINT(bool        , USETOOCOARSESAMPLING, 1         , UseTooCoarseSampling         ) SEP \
    ELTINT(int         , WIDTHMASKEDGE      , 0         , WidthMaskEdge                ) SEP \
    ELTINT(bool        , ISFLIP             , 1         , IsFlip                       ) SEP \
    ELTINT(int         , ORIENTATIONSID     , 0         , OrientationsID               ) SEP \
    ELTDOU(FDOUBLE     , ANGLEPSIFLIPRATIO  , 0.0       , AnglePsiFlipRatio            ) SEP \
    ELTDOU(FDOUBLE     , AUTOPICKFIGUREOFMERIT, 0.0       , AutopickFigureOfMerit        ) SEP \
    ELTINT(int         , HELICALTUBEID      , 0         , HelicalTubeID                ) SEP \
    ELTDOU(FDOUBLE     , HELICALTUBEPITCH   , 0.0       , HelicalTubePitch             ) SEP \
    ELTDOU(FDOUBLE     , HELICALTRACKLENGTH , 0.0       , HelicalTrackLength           ) SEP \
    ELTINT(int         , PARTICLEID         , 0         , ParticleId                   ) SEP \
    ELTDOU(FDOUBLE     , PARTICLEFIGUREOFMERIT, 0.0       , ParticleFigureOfMerit        ) SEP \
    ELTDOU(FDOUBLE     , KULLBACKLEIBLERDIVERGENCE, 0.0       , KullbackLeiblerDivergence    ) SEP \
    ELTINT(int         , BEAMTILTCLASS      , 0         , BeamTiltClass                ) SEP \
    ELTSTR(std::string , PARTICLENAME       , ""        , ParticleName                 ) SEP \
    ELTSTR(std::string , ORIGINALPARTICLENAME, ""        , OriginalParticleName         ) SEP \
    ELTINT(int         , NROFFRAMES         , 0         , NrOfFrames                   ) SEP \
    ELTINT(int         , AVERAGENROFFRAMES  , 0         , AverageNrOfFrames            ) SEP \
    ELTINT(int         , MOVIEFRAMESRUNNINGAVERAGE, 0         , MovieFramesRunningAverage    ) SEP \
    ELTINT(int         , PARTICLENUMBER     , 0         , ParticleNumber               ) SEP \
    ELTINT(int         , PIPELINEJOBCOUNTER , 0         , PipeLineJobCounter           ) SEP \
    ELTSTR(std::string , PIPELINENODENAME   , ""        , PipeLineNodeName             ) SEP \
    ELTINT(int         , PIPELINENODETYPE   , 0         , PipeLineNodeType             ) SEP \
    ELTSTR(std::string , PIPELINEPROCESSALIAS     , ""        ,  PipeLineProcessAlias        ) SEP \
    ELTSTR(std::string , PIPELINEPROCESSNAME       , ""        ,  PipeLineProcessName         ) SEP \
    ELTINT(int         , PIPELINEPROCESSTYPE, 0         , PipeLineProcessType          ) SEP \
    ELTINT(int         , PIPELINEPROCESSSTATUS, 0         , PipeLineProcessStatus        ) SEP \
    ELTSTR(std::string , PIPELINEEDGEFROMNODE , ""        , PipeLineEdgeFromNode                         ) SEP \
    ELTSTR(std::string , PIPELINEEDGETONODE , ""        , PipeLineEdgeToNode           ) SEP \
    ELTSTR(std::string , PIPELINEEDGEPROCESS, ""        , PipeLineEdgeProcess          ) SEP \
    ELTDOU(FDOUBLE     , FINALRESOLUTION    , 0.0       , FinalResolution              ) SEP \
    ELTDOU(FDOUBLE     , BFACTORUSEDFORSHARPENING, 0.0       , BfactorUsedForSharpening     ) SEP \
    ELTDOU(FDOUBLE     , FOURIERSHELLCORRELATION, 0.0       , FourierShellCorrelation      ) SEP \
    ELTDOU(FDOUBLE     , FOURIERSHELLCORRELATIONCORRECTED, 0.0       , FourierShellCorrelationCorrected ) SEP \
    ELTDOU(FDOUBLE     , FOURIERSHELLCORRELATIONMASKEDMAPS, 0.0       , FourierShellCorrelationMaskedMaps ) SEP \
    ELTDOU(FDOUBLE     , FOURIERSHELLCORRELATIONUNMASKEDMAPS, 0.0       , FourierShellCorrelationUnmaskedMaps ) SEP \
	/* bug here */ \
	ELTDOU(FDOUBLE     , CORRECTEDFOURIERSHELLCORRELATIONPHASERANDOMIZEDMASKEDMAPS, 0.0       , CorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps ) SEP \
    ELTDOU(FDOUBLE     , AMPLITUDECORRELATIONMASKEDMAPS, 0.0       , AmplitudeCorrelationMaskedMaps ) SEP \
    ELTDOU(FDOUBLE     , AMPLITUDECORRELATIONUNMASKEDMAPS, 0.0       , AmplitudeCorrelationUnmaskedMaps ) SEP \
    ELTDOU(FDOUBLE     , DIFFERENTIALPHASERESIDUALMASKEDMAPS, 0.0       , DifferentialPhaseResidualMaskedMaps ) SEP \
    ELTDOU(FDOUBLE     , DIFFERENTIALPHASERESIDUALUNMASKEDMAPS, 0.0       , DifferentialPhaseResidualUnmaskedMaps ) SEP \
    ELTDOU(FDOUBLE     , FITTEDINTERCEPTGUINIERPLOT, 0.0       , FittedInterceptGuinierPlot   ) SEP \
    ELTDOU(FDOUBLE     , FITTEDSLOPEGUINIERPLOT, 0.0       , FittedSlopeGuinierPlot       ) SEP \
    ELTDOU(FDOUBLE     , CORRELATIONFITGUINIERPLOT, 0.0       , CorrelationFitGuinierPlot    ) SEP \
    ELTDOU(FDOUBLE     , LOGAMPLITUDESORIGINAL, 0.0       , LogAmplitudesOriginal        ) SEP \
    ELTDOU(FDOUBLE     , LOGAMPLITUDESMTFCORRECTED, 0.0       , LogAmplitudesMTFCorrected    ) SEP \
    ELTDOU(FDOUBLE     , LOGAMPLITUDESWEIGHTED, 0.0       , LogAmplitudesWeighted        ) SEP \
    ELTDOU(FDOUBLE     , LOGAMPLITUDESSHARPENED, 0.0       , LogAmplitudesSharpened       ) SEP \
    ELTDOU(FDOUBLE     , LOGAMPLITUDESINTERCEPT, 0.0       , LogAmplitudesIntercept       ) SEP \
    ELTDOU(FDOUBLE     , RESOLUTIONSQUARED  , 0.0       , ResolutionSquared            ) SEP \
    ELTDOU(FDOUBLE     , MTFVALUE           , 0.0       , MtfValue                     ) SEP \
    ELTDOU(FDOUBLE     , RANDOMISEFROM      , 0.0       , RandomiseFrom                ) SEP \
    ELTSTR(std::string , UNFILTEREDMAPHALF1 , ""        , UnfilteredMapHalf1           ) SEP \
    ELTSTR(std::string , UNFILTEREDMAPHALF2 , ""        , UnfilteredMapHalf2           ) SEP \
    ELTINT(bool        , IS3DSAMPLING       , 1         , Is3DSampling                 ) SEP \
    ELTINT(bool        , IS3DTRANSLATIONALSAMPLING, 1         , Is3DTranslationalSampling    ) SEP \
    ELTINT(int         , HEALPIXORDER       , 0         , HealpixOrder                 ) SEP \
    ELTDOU(FDOUBLE     , TILTANGLELIMIT     , 0.0       , TiltAngleLimit               ) SEP \
    ELTDOU(FDOUBLE     , OFFSETRANGE        , 0.0       , OffsetRange                  ) SEP \
    ELTDOU(FDOUBLE     , OFFSETSTEP         , 0.0       , OffsetStep                   ) SEP \
    ELTDOU(FDOUBLE     , HELICALOFFSETSTEP  , 0.0       , HelicalOffsetStep            ) SEP \
    ELTDOU(FDOUBLE     , SAMPLINGPERTURBINSTANCE, 0.0       , SamplingPerturbInstance      ) SEP \
    ELTDOU(FDOUBLE     , SAMPLINGPERTURBFACTOR, 0.0       , SamplingPerturbFactor        ) SEP \
    ELTDOU(FDOUBLE     , PSISTEP            , 0.0       , PsiStep                      ) SEP \
    ELTSTR(std::string , SYMMETRYGROUP      , ""        , SymmetryGroup                ) SEP \
    ELTINT(int         , SELECTED           , 0         , Selected                     ) SEP \
    ELTDOU(FDOUBLE     , PARTICLESELECTZSCORE, 0.0       , ParticleSelectZScore         ) SEP \
    ELTINT(int         , SORTEDINDEX        , 0         , SortedIndex                  ) SEP \
    ELTSTR(std::string , STARFILEMOVIEPARTICLES, ""        , StarFileMovieParticles       ) SEP \
    ELTDOU(FDOUBLE     , PERFRAMECUMULATIVEWEIGHT, 0.0       , PerFrameCumulativeWeight     ) SEP \
    ELTDOU(FDOUBLE     , PERFRAMERELATIVEWEIGHT, 0.0       , PerFrameRelativeWeight       ) SEP \
    ELTDOU(FDOUBLE     , RESOLUTION         , 0.0       , Resolution                   ) SEP \
    ELTDOU(FDOUBLE     , ANGSTROMRESOLUTION , 0.0       , AngstromResolution           ) SEP \
    ELTDOU(FDOUBLE     , RESOLUTIONINVERSEPIXEL, 0.0       , ResolutionInversePixel       ) SEP \
    ELTINT(int         , SPECTRALINDEX      , 0         , SpectralIndex                )



enum MetaDataLabel{
#define SEP ,
#define ELTINT(T,N,I,C) C
#define ELTDOU(T,N,I,C) C
#define ELTSTR(T,N,I,C) C
    METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    ,UnDefinedLabel
};

namespace LabelParser{
    MetaDataLabel stringToLabel(const std::string& name);
    std::string labelToString(MetaDataLabel label);
};

// one metadata element
struct MetaDataElem{
    size_t INNERID = -1;
#define SEP
#define ELTINT(T,N,I,C) T N = I;
#define ELTDOU(T,N,I,C) T N = I;
#define ELTSTR(T,N,I,C) T N = I;
    METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    
    void set(MetaDataLabel label,std::string val){
        switch(label)
        {
#define SEP
#define ELTINT(T,N,I,C) case C : N = atoi(val.c_str()); break;
#define ELTDOU(T,N,I,C) case C : N = atof(val.c_str()); break;
#define ELTSTR(T,N,I,C) case C : N = val            ; break;
            METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
        }
    }
    
    void get(MetaDataLabel label,std::ostream& os) const {
        // os.precision(6);
        auto printDouble = [&](double v){
            if((fabs(v) > 0. && fabs(v) < 0.001) || fabs(v) > 100000.)
                os<<std::setw(12)<<std::scientific<<v<<"\t";
            else
                os<<std::setw(12)<<std::fixed<<v<<"\t";
        };
        auto printInt = [&](int v){
            os<<std::setw(6)<<v<<"\t";
        };
        
        switch(label)
        {
#define SEP
#define ELTINT(T,N,I,C) case C : printInt(N)        ; break;
#define ELTDOU(T,N,I,C) case C : printDouble(N)     ; break;
#define ELTSTR(T,N,I,C) case C : os<<N<<"\t";  ; break;
            METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
        }
    }
    
    void operator=(MetaDataElem const & rhs) {
        INNERID = rhs.INNERID;
#define SEP
#define ELTINT(T,N,I,C) N = rhs.N;
#define ELTDOU(T,N,I,C) N = rhs.N;
#define ELTSTR(T,N,I,C) N = rhs.N;
        METADATA_ELTS
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    }

    bool operator==(MetaDataElem const & rhs) const {
        const char* diff = NULL;
		//
		// NR_SIGN is a challenge, because it is an integer but the exact number (number of significant samples) is imprecise
		//
		auto compareInts = [&](const char* name, int lhs, int rhs) {
			if (lhs == rhs) return;
			if (strcmp(name, "NR_SIGN") == 0 && std::abs(lhs-rhs) < std::max(lhs,rhs)/100) return;
			diff = name;
			std::cerr << "MetaDataElem:" << name << " lhs:" << lhs << " != " << rhs<< std::endl;
		};
#define SEP
#define ELTINT(T,N,I,C) compareInts(#N, N, rhs.N);
#define ELTDOU(T,N,I,C) if (!nearEnoughTemplate(N, rhs.N)) { diff = #N;	std::cerr << "MetaDataElem:" << #N << " lhs:" << N << " != " << rhs.N << std::endl; }
#define ELTSTR(T,N,I,C) if (!nearEnoughTemplate(N, rhs.N)) { diff = #N;	std::cerr << "MetaDataElem:" << #N << " lhs:" << N << " != " << rhs.N << std::endl; }
    METADATA_ELTS
        if (INNERID != rhs.INNERID) { diff = "innerID"; std::cerr << "MetaDataElem:INNERID lhs:" << INNERID << " != " << rhs.INNERID << std::endl; }
        if (!diff) {
			return true;
		}
        return false;
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    }

    void print(std::ostream& os) const {
        os << "MetaDataElems{"
#define SEP
#define ELTINT(T,N,I,C) << " " << #N << ":" << N
#define ELTDOU(T,N,I,C) << " " << #N << ":" << N
#define ELTSTR(T,N,I,C) << " " << #N << ":" << N
    METADATA_ELTS
        << " }";
#undef ELTINT
#undef ELTDOU
#undef ELTSTR
#undef SEP
    }
    
};

inline std::ostream& operator<<(std::ostream& os, MetaDataElem const & rhs) { rhs.print(os); return os; }

class MetaDataTable : public NoCopy {
public:
    MetaDataTable():N(0),latestLabelStatus(nullptr) {}
    ~MetaDataTable() { freeMetaData(); }
    void clear()     { freeMetaData(); }
    
    // initialize the metadata table
    void readFromMetaDataElements(std::vector<MetaDataElem>& _metaDataElems);

	// Read/Write *.star file from/to metadata
    void readFromStar(std::string fn_star);
    void writeToStar (std::string fn_star) const;

	//
    void printTable(std::ostream& os = std::cout) const;
    
    // get the metadata elements size
    int numberOfParticles(int random_subset = -1) const {
        if (random_subset==1) return nr_ori_particles_subset1;
        else if (random_subset==2) return nr_ori_particles_subset2;
        else return N;
    }
    void selectHalfMetadata(int random_subset) {
        if (random_subset==1) {
            subset_start = 0;subset_end = nr_ori_particles_subset1;
        }
        else if (random_subset==2) {
            subset_start = nr_ori_particles_subset1;subset_end = nr_ori_particles_subset1+nr_ori_particles_subset2;
        }
        else {
            subset_start = 0;subset_end = N;
        }
    }
    int numberOfGroups()		const 	{return nr_groups;}
    int numberOfMicrographs()	const	{return nr_micrograph;}
    //
    void append(std::vector<MetaDataElem>& _metadataElems);
    
    inline MetaDataElem       & operator[](int i)       {assert(i>-1);assert(i<(subset_end-subset_start)); return metaDataElems[metadata_order[i+subset_start]];}
    inline MetaDataElem const & operator[](int i) const {assert(i>-1);assert(i<(subset_end-subset_start)); return metaDataElems[metadata_order[i+subset_start]];}
    inline MetaDataElem 	  & accessAll(int i)		{assert(i>-1);assert(i<N); return metaDataElems[metadata_order[i]];}
    inline MetaDataElem const & accessAll(int i) const  {assert(i>-1);assert(i<N); return metaDataElems[metadata_order[i]];}
    // randomly shuffle the metadata elements
    // and split metadata to two random subset if no randomsubset read from *.star(only used in autorefine)
    void shuffle(int random_seed = -1,bool do_split_random_halves = false);
    
    // some statistics function
    void statClassImageNumber(std::vector<int> &imageNumberPerClass,bool showHist = false);
    // some fliter function
    void fliterByClass(int selectedClass);
    //
    bool containLabel(MetaDataLabel label);
    
    void unitTest(std::string fn);
    
private:
    void freeMetaData();

	// all
	class LabelStatus;
	LabelStatus*			  latestLabelStatus;
    std::vector<MetaDataElem> metaDataElems;

    // other metadata elements
    std::map<int,std::string>              undefinedMetaDataLabels; // label:column_index
    std::map<int,std::vector<std::string>> undefinedMetaDataElems;	// inner_id:....

    // metadata element number
    int N;
    int nr_groups		= 0;
    int nr_micrograph	= 0;
    // bool do_split_random_halves = false; // split data to two random subset(in autorefine)
    int nr_ori_particles_subset1 = 0;
    int nr_ori_particles_subset2 = 0;
    int subset_start,subset_end;
	//
    std::vector<int> metadata_order;

	int getImageNumber(std::string starFileName);

    void printTable(std::ostream& os, LabelStatus const& labelStatus) const;

};

enum ElemType {ElemTypeChar,ElemTypeInt,ElemTypeDouble};

class SmallMetataDataTable : public NoCopy {
private:
    std::string tableName;
    std::vector<std::string> MetaDataElemsName;
    std::vector<ElemType> MetaDataElemsType;
    std::vector<void*> MetaDataElems;
    int MetaDataElemsNumber;
public:
    SmallMetataDataTable(std::string _tableName){tableName = _tableName;}
    ~SmallMetataDataTable(){
        MetaDataElemsName.resize(0);
        MetaDataElemsType.resize(0);
        MetaDataElems.resize(0);
        MetaDataElemsNumber = 0;
    }
    void appendName(std::initializer_list<std::string> _MetaDataElemsName){
        for (auto const &__name : _MetaDataElemsName)
            MetaDataElemsName.push_back(__name);
    }
    void appendType(std::initializer_list<ElemType> _MetaDataElemsType){
        for (auto const &__type : _MetaDataElemsType)
            MetaDataElemsType.push_back(__type);
    }
    void appendElem(std::initializer_list<void*> _MetaDataElems,int _MetaDataElemsNumber){
        MetaDataElemsNumber = _MetaDataElemsNumber;
        for (auto const &__elem : _MetaDataElems)
            MetaDataElems.push_back(__elem);
    }
    void print(std::ostream& os)
    {
        assert(MetaDataElemsName.size()==MetaDataElemsType.size());
        assert(MetaDataElemsType.size()==MetaDataElems.size());
        os << tableName <<std::endl<<std::endl<<"loop_"<<std::endl;
		int ScaleCorrectionIndex = 0; // remove rlnScaleCorrection
        for (int i = 0; i < MetaDataElemsName.size(); i++) {
			if (MetaDataElemsName[i].compare("ScaleCorrection") == 0 ) { // remove rlnScaleCorrection
				ScaleCorrectionIndex = i;
				continue;
			}
			else if (i > ScaleCorrectionIndex)
				os << "_rln" << MetaDataElemsName[i] << " #" << std::to_string((long long)i) << std::endl;
            else
				os << "_rln"<<MetaDataElemsName[i]<<" #"<<std::to_string((long long)i+1)<<std::endl;
        }
        for (int metadataIndex = 0; metadataIndex < MetaDataElemsNumber; metadataIndex++) {
            for (int i = 0; i < MetaDataElems.size(); i++) {
				if (i == ScaleCorrectionIndex)  continue;
                switch (MetaDataElemsType[i]) {
                    case ElemTypeChar: {
                        auto string_ptr = (std::string**)MetaDataElems[i];
                        os << std::setw(12) << *string_ptr[metadataIndex] << " ";
                    }   break;
                    case ElemTypeInt: {
                        auto int_ptr = (int*)MetaDataElems[i];
                        os << std::setw(12) << int_ptr[metadataIndex] << " ";
                    }   break;
                    case ElemTypeDouble: {
                        auto* double_ptr = (double*)MetaDataElems[i];
                        auto var = double_ptr[metadataIndex];
                        if (var < 1 && var != 0) os << std::setw(12) << std::scientific << std::setprecision(7) << var << " ";
                        else os << std::setw(12) << std::fixed << std::setprecision(7) << var << " ";
                    }   break;
                    default:
                        break;
                }
            }
            os << std::endl;
        }
        os << std::endl << std::endl;
    }
    bool read(std::istream& is)
    {
        assert(MetaDataElemsName.size()==MetaDataElemsType.size());
        assert(MetaDataElemsType.size()==MetaDataElems.size());
        // reset to beginning
        is.clear();is.seekg(0, std::ios::beg);
        bool startingRead = false;
        std::string line;
        int metadataIndex = 0;
        //
        while (true)
        {
            getline(is,line);
            if ( (line.find(tableName) !=std::string::npos) ){
                startingRead = true;
                getline(is,line);assert(line=="");// escape a empty line
                getline(is,line);assert(line.find("loop_")!=std::string::npos);// escape loop_ line
                // head
                for (int headIndex = 0; headIndex < MetaDataElems.size(); headIndex++) {
                    getline(is,line);
                    if( line.find(MetaDataElemsName[headIndex]) == std::string::npos )
                        std::cerr<<"cannot find "+MetaDataElemsName[headIndex]+" in _model.star."<<std::endl;
                }
                break;
            }
            ERROR_CHECK(is.eof(), "end of model file,can not find "+tableName+".");
        }
        //
        bool metadata_num_fit = true;
        for (int metadataIndex = 0; metadataIndex < MetaDataElemsNumber; metadataIndex++)
        {
            getline(is,line);
            if (line.empty()) {
                metadata_num_fit = false;
                break;
            }
            std::istringstream lineStream(line);
            for (int i = 0; i < MetaDataElems.size(); i++)
            {
                switch (MetaDataElemsType[i]) {
                    case ElemTypeChar:{
                        std::string lineTmp;lineStream >> lineTmp;
                        ((std::string**)MetaDataElems[i])[metadataIndex] = 
#include "./util_heap_undefs.h"
							new std::string(lineTmp);
#include "./util_heap_defs.h"
					}	break;
                    case ElemTypeInt:
                        lineStream >> ((int*)MetaDataElems[i])[metadataIndex];
                    	break;
                    case ElemTypeDouble:
                    	lineStream >> ((double*)MetaDataElems[i])[metadataIndex];
                        break;
                    default:
                        break;
                }
            }
        }
        //if (!metadata_num_fit) {MASTERNODE std::cout<<"_rln"<<MetaDataElemsName[0]<<" has less elements to read!!!!!"<<std::endl;}
        return metadata_num_fit;
    }
    bool contain(std::istream& is,std::string head)
    {
        assert(MetaDataElemsName.size()==MetaDataElemsType.size());
        assert(MetaDataElemsType.size()==MetaDataElems.size());
        // reset to beginning
        is.clear();is.seekg(0, std::ios::beg);
        bool startingRead = false;
        std::string line;
        bool contain = false;
        //
        while (true)
        {
            getline(is,line);
            if ( (line.find(tableName) != std::string::npos) ){
                startingRead = true;
                getline(is,line);assert(line=="");// escape a empty line
                getline(is,line);assert(line.find("loop_")!=std::string::npos);// escape loop_ line
                // head
                while (startingRead) {
                    getline(is,line);
                    if(line.find("_rln") == std::string::npos) break;
                    if(line.find(head) != std::string::npos) contain = true;
                }
                break;
            }
            ERROR_CHECK(is.eof(), "end of model file,can not find "+tableName+".");
        }
        return contain;
    }
};

#endif
