# Makefile created by mkmf $Id: mkmf,v 16.1 2010/05/19 18:49:19 fms Exp $ 


include ../../Rules.make

include mkmf.template.ifc


.DEFAULT:
	-echo $@ does not exist.
all: libtopaz.a
FMS_coupler_util.o: ./FMS_coupler_util.F90 mpp.o coupler_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./FMS_coupler_util.F90
FMS_ocmip2_co2calc.o: ./FMS_ocmip2_co2calc.F90
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./FMS_ocmip2_co2calc.F90
algebraiclength.o: ./algebraiclength.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./algebraiclength.F90
alpha_mnb.o: ./alpha_mnb.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./alpha_mnb.F90
atmos_ocean_fluxes.o: ./atmos_ocean_fluxes.F90 mpp.o coupler_types.o field_manager.o fm_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./atmos_ocean_fluxes.F90
axis_utils.o: ./axis_utils.F90 mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./axis_utils.F90
cmue_c.o: ./cmue_c.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./cmue_c.F90
cmue_d.o: ./cmue_d.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./cmue_d.F90
cmue_ma.o: ./cmue_ma.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./cmue_ma.F90
cmue_rf.o: ./cmue_rf.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./cmue_rf.F90
cmue_sg.o: ./cmue_sg.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./cmue_sg.F90
compute_cpsi3.o: ./compute_cpsi3.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./compute_cpsi3.F90
compute_rist.o: ./compute_rist.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./compute_rist.F90
constants.o: ./constants.F90
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./constants.F90
coupler_types.o: ./coupler_types.F90 field_manager.o mpp.o fm_util.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./coupler_types.F90
create_xgrid.o: ./create_xgrid.c
	$(CC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	./create_xgrid.c
data_override.o: ./data_override.F90 platform.o constants.o mpp_io.o mpp.o horiz_interp.o time_interp_external.o fms_io.o fms.o axis_utils.o mpp_domains.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./data_override.F90
diag_axis.o: ./diag_axis.F90 mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./diag_axis.F90
diag_data.o: ./diag_data.F90 time_manager.o mpp_domains.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./diag_data.F90
diag_grid.o: ./diag_grid.F90 constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./diag_grid.F90
diag_manager.o: ./diag_manager.F90 time_manager.o mpp_io.o mpp.o fms.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o constants.o mpp_domains.o fms_io.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./diag_manager.F90
diag_output.o: ./diag_output.F90 mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./diag_output.F90
diag_table.o: ./diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./diag_table.F90
diag_util.o: ./diag_util.F90 diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./diag_util.F90
diff_face.o: ./diff_face.F90 util.o tridiagonal_gotm.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./diff_face.F90
dissipationeq.o: ./dissipationeq.F90 turbulence.o util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./dissipationeq.F90
epsbalgebraic.o: ./epsbalgebraic.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./epsbalgebraic.F90
field_manager.o: ./field_manager.F90 mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./field_manager.F90
fk_craig.o: ./fk_craig.F90
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./fk_craig.F90
fm_util.o: ./fm_util.F90 field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./fm_util.F90
fms.o: ./fms.F90 mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./fms.F90
fms_io.o: ./fms_io.F90 mpp_io.o mpp_domains.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./fms_io.F90
generic_BLING.o: ./generic_BLING.F90 coupler_types.o field_manager.o mpp.o fms.o time_manager.o fm_util.o diag_manager.o constants.o generic_tracer_utils.o FMS_ocmip2_co2calc.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./generic_BLING.F90
generic_CFC.o: ./generic_CFC.F90 coupler_types.o field_manager.o mpp.o time_manager.o fm_util.o generic_tracer_utils.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./generic_CFC.F90
generic_ERGOM.o: ./generic_ERGOM.F90 coupler_types.o field_manager.o fms.o mpp.o time_manager.o fm_util.o diag_manager.o constants.o generic_tracer_utils.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./generic_ERGOM.F90
generic_TOPAZ.o: ./generic_TOPAZ.F90 coupler_types.o field_manager.o mpp.o time_manager.o fm_util.o diag_manager.o constants.o generic_tracer_utils.o FMS_ocmip2_co2calc.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./generic_TOPAZ.F90
generic_TOPAZ_model.o: ./generic_TOPAZ_model.F90 ocean_types.o generic_tracer_utils.o time_manager.o coupler_types.o netcdf_var.o netcdf_io.o fms_io.o ocean_barotropic.o ocean_grids.o ocean_domains.o ocean_thickness.o ocean_workspace.o ocean_topog.o ocean_util.o ocean_obc.o ocean_operators.o ocean_coriolis.o ocean_velocity.o diag_manager.o ocean_tracer.o ocean_tracer_util.o ocean_generic_tracer.o generic_tracer.o generic_TOPAZ.o time_interp.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./generic_TOPAZ_model.F90
generic_tracer.o: ./generic_tracer.F90 fms.o field_manager.o mpp.o time_manager.o coupler_types.o generic_tracer_utils.o generic_CFC.o generic_ERGOM.o generic_TOPAZ.o generic_BLING.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./generic_tracer.F90
generic_tracer_utils.o: ./generic_tracer_utils.F90 coupler_types.o FMS_coupler_util.o atmos_ocean_fluxes.o mpp.o time_manager.o diag_manager.o field_manager.o fms.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./generic_tracer_utils.F90
genericeq.o: ./genericeq.F90 turbulence.o util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./genericeq.F90
get_cal_time.o: ./get_cal_time.F90 fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./get_cal_time.F90
horiz_interp.o: ./horiz_interp.F90 fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o mpp_io.o mpp_domains.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./horiz_interp.F90
horiz_interp_bicubic.o: ./horiz_interp_bicubic.F90 mpp.o fms.o horiz_interp_type.o constants.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./horiz_interp_bicubic.F90
horiz_interp_bilinear.o: ./horiz_interp_bilinear.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./horiz_interp_bilinear.F90
horiz_interp_conserve.o: ./horiz_interp_conserve.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./horiz_interp_conserve.F90
horiz_interp_spherical.o: ./horiz_interp_spherical.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./horiz_interp_spherical.F90
horiz_interp_type.o: ./horiz_interp_type.F90 mpp.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./horiz_interp_type.F90
internal_wave.o: ./internal_wave.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./internal_wave.F90
ispralength.o: ./ispralength.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ispralength.F90
kbalgebraic.o: ./kbalgebraic.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./kbalgebraic.F90
kbeq.o: ./kbeq.F90 turbulence.o util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./kbeq.F90
lengthscaleeq.o: ./lengthscaleeq.F90 turbulence.o util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./lengthscaleeq.F90
memutils.o: ./memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./memutils.F90
mosaic.o: ./mosaic.F90 mpp.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./mosaic.F90
mosaic_util.o: ./mosaic_util.c
	$(CC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	./mosaic_util.c
mpp.o: ./mpp.F90 mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./mpp.F90
mpp_data.o: ./mpp_data.F90 mpp_parameter.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./mpp_data.F90
mpp_domains.o: ./mpp_domains.F90 mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./mpp_domains.F90
mpp_io.o: ./mpp_io.F90 mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./mpp_io.F90
mpp_memutils.o: ./mpp_memutils.F90 mpp.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./mpp_memutils.F90
mpp_parameter.o: ./mpp_parameter.F90
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./mpp_parameter.F90
mpp_pset.o: ./mpp_pset.F90 mpp.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./mpp_pset.F90
netcdf_io.o: ./netcdf_io.F90 netcdf_var.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./netcdf_io.F90
netcdf_var.o: ./netcdf_var.F90
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./netcdf_var.F90
ocean_age_tracer.o: ./ocean_age_tracer.F90 field_manager.o mpp.o time_manager.o ocean_tpm_util.o fm_util.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_age_tracer.F90
ocean_barotropic.o: ./ocean_barotropic.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o time_manager.o time_interp_external.o ocean_bih_friction.o ocean_domains.o ocean_lap_friction.o ocean_obc_barotrop.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_barotropic.F90
ocean_bgc_restore.o: ./ocean_bgc_restore.F90 time_manager.o diag_manager.o field_manager.o fms.o fms_io.o mpp.o time_interp_external.o mpp_domains.o constants.o ocean_tpm_util.o fm_util.o coupler_types.o ocean_types.o ocmip2_co2calc.o ocean_util.o atmos_ocean_fluxes.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_bgc_restore.F90
ocean_bih_friction.o: ./ocean_bih_friction.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_bihcgrid_friction.o ocean_bihcst_friction.o ocean_bihgen_friction.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_bih_friction.F90
ocean_bih_tracer.o: ./ocean_bih_tracer.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_operators.o ocean_obc.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_bih_tracer.F90
ocean_bihcgrid_friction.o: ./ocean_bihcgrid_friction.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_bihcgrid_friction.F90
ocean_bihcst_friction.o: ./ocean_bihcst_friction.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_bihcst_friction.F90
ocean_bihgen_friction.o: ./ocean_bihgen_friction.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_bihgen_friction.F90
ocean_blob.o: ./ocean_blob.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_blob_diag.o ocean_blob_dynamic_bottom.o ocean_blob_dynamic_free.o ocean_blob_static_free.o ocean_blob_static_bottom.o ocean_blob_util.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_blob.F90
ocean_blob_diag.o: ./ocean_blob_diag.F90 constants.o fms.o mpp.o ocean_blob_util.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_blob_diag.F90
ocean_blob_dynamic_bottom.o: ./ocean_blob_dynamic_bottom.F90 constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_blob_util.o ocean_density.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_blob_dynamic_bottom.F90
ocean_blob_dynamic_free.o: ./ocean_blob_dynamic_free.F90 constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_blob_util.o ocean_density.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_blob_dynamic_free.F90
ocean_blob_static_bottom.o: ./ocean_blob_static_bottom.F90 fms.o mpp_domains.o mpp.o ocean_blob_util.o ocean_density.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_blob_static_bottom.F90
ocean_blob_static_free.o: ./ocean_blob_static_free.F90 fms.o ocean_blob_util.o ocean_density.o ocean_types.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_blob_static_free.F90
ocean_blob_util.o: ./ocean_blob_util.F90 constants.o fms.o mpp_domains.o mpp.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_blob_util.F90
ocean_convect.o: ./ocean_convect.F90 constants.o diag_manager.o fms.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_convect.F90
ocean_coriolis.o: ./ocean_coriolis.F90 constants.o diag_manager.o fms.o mpp.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_coriolis.F90
ocean_density.o: ./ocean_density.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o platform.o time_manager.o field_manager.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_pressure.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_density.F90
ocean_domains.o: ./ocean_domains.F90 fms.o mpp_domains.o mpp.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_domains.F90
ocean_form_drag.o: ./ocean_form_drag.F90 constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_form_drag.F90
ocean_frazil.o: ./ocean_frazil.F90 constants.o diag_manager.o fms.o mpp.o mpp_domains.o ocean_domains.o ocean_parameters.o ocean_tpm_util.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_frazil.F90
ocean_generic_tracer.o: ./ocean_generic_tracer.F90 ocean_tpm_util.o ocean_types.o ocean_parameters.o mpp.o field_manager.o generic_tracer.o generic_tracer_utils.o coupler_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_generic_tracer.F90
ocean_grids.o: ./ocean_grids.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o mosaic.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_grids.F90
ocean_ibgc.o: ./ocean_ibgc.F90 diag_manager.o field_manager.o mpp.o fms.o fms_io.o time_manager.o time_interp_external.o mpp_domains.o constants.o ocean_tpm_util.o fm_util.o ocean_types.o ocmip2_co2calc.o coupler_types.o atmos_ocean_fluxes.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_ibgc.F90
ocean_lap_friction.o: ./ocean_lap_friction.F90 diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_domains.o ocean_lapcst_friction.o ocean_lapgen_friction.o ocean_lapcgrid_friction.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_lap_friction.F90
ocean_lap_tracer.o: ./ocean_lap_tracer.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_lap_tracer.F90
ocean_lapcgrid_friction.o: ./ocean_lapcgrid_friction.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_lapcgrid_friction.F90
ocean_lapcst_friction.o: ./ocean_lapcst_friction.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_lapcst_friction.F90
ocean_lapgen_friction.o: ./ocean_lapgen_friction.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_lapgen_friction.F90
ocean_momentum_source.o: ./ocean_momentum_source.F90 constants.o diag_manager.o field_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_momentum_source.F90
ocean_obc.o: ./ocean_obc.F90 constants.o data_override.o diag_manager.o fms.o fms_io.o mpp_io.o mpp_domains.o mpp.o time_interp_external.o time_manager.o tracer_manager.o ocean_util.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_obc_barotrop.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_obc.F90
ocean_obc_barotrop.o: ./ocean_obc_barotrop.F90 constants.o data_override.o diag_manager.o fms.o fms_io.o mpp_io.o mpp_domains.o mpp.o time_interp_external.o time_manager.o tracer_manager.o ocean_util.o ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_obc_barotrop.F90
ocean_operators.o: ./ocean_operators.F90 mpp_domains.o mpp.o fms.o ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_operators.F90
ocean_parameters.o: ./ocean_parameters.F90 mpp.o fms.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_parameters.F90
ocean_passive.o: ./ocean_passive.F90 constants.o field_manager.o fm_util.o fms.o mpp.o ocean_domains.o ocean_tpm_util.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_passive.F90
ocean_pert_co2.o: ./ocean_pert_co2.F90 atmos_ocean_fluxes.o diag_manager.o field_manager.o time_manager.o fms.o mpp.o mpp_domains.o constants.o ocean_tpm_util.o fm_util.o coupler_types.o ocean_types.o ocmip2_co2calc.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_pert_co2.F90
ocean_po4_pre.o: ./ocean_po4_pre.F90 diag_manager.o field_manager.o mpp.o fms.o time_manager.o mpp_domains.o ocean_tpm_util.o fm_util.o ocean_types.o ocean_tracer_diag.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_po4_pre.F90
ocean_pressure.o: ./ocean_pressure.F90 constants.o diag_manager.o fms.o mpp_io.o mpp.o mpp_domains.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o ocean_obc.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_pressure.F90
ocean_residency.o: ./ocean_residency.F90 time_manager.o field_manager.o mpp.o diag_manager.o ocean_tpm_util.o fm_util.o ocean_types.o ocean_residency_meta.o ocean_residency_ml.o ocean_residency_range.o ocean_residency_integrand.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_residency.F90
ocean_residency_integrand.o: ./ocean_residency_integrand.F90 field_manager.o mpp.o fms.o ocean_residency_meta.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_residency_integrand.F90
ocean_residency_meta.o: ./ocean_residency_meta.F90 field_manager.o fm_util.o mpp.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_residency_meta.F90
ocean_residency_ml.o: ./ocean_residency_ml.F90 field_manager.o mpp.o ocean_residency_meta.o ocean_types.o ocean_tracer_diag.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_residency_ml.F90
ocean_residency_range.o: ./ocean_residency_range.F90 field_manager.o mpp.o ocean_residency_meta.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_residency_range.F90
ocean_shortwave.o: ./ocean_shortwave.F90 constants.o diag_manager.o fms.o mpp.o ocean_domains.o ocean_parameters.o ocean_shortwave_csiro.o ocean_shortwave_jerlov.o ocean_shortwave_gfdl.o ocean_types.o ocean_tpm_util.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_shortwave.F90
ocean_shortwave_csiro.o: ./ocean_shortwave_csiro.F90 axis_utils.o field_manager.o fms.o mpp.o time_interp_external.o ocean_domains.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_shortwave_csiro.F90
ocean_shortwave_gfdl.o: ./ocean_shortwave_gfdl.F90 axis_utils.o constants.o diag_manager.o field_manager.o fms.o mpp.o time_interp_external.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_shortwave_gfdl.F90
ocean_shortwave_jerlov.o: ./ocean_shortwave_jerlov.F90 constants.o fms.o mpp.o time_interp_external.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_shortwave_jerlov.F90
ocean_tempsalt.o: ./ocean_tempsalt.F90 constants.o fms.o mpp.o ocean_domains.o ocean_parameters.o ocean_tpm_util.o ocean_types.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_tempsalt.F90
ocean_thickness.o: ./ocean_thickness.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_domains.o ocean_grids.o ocean_parameters.o ocean_tracer_util.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_thickness.F90
ocean_topog.o: ./ocean_topog.F90 fms.o mpp_domains.o mpp.o axis_utils.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_topog.F90
ocean_tpm.o: ./ocean_tpm.F90 mpp.o mpp_domains.o ocean_types.o ocean_age_tracer.o ocean_residency.o ocean_pert_co2.o ocmip2_abiotic.o ocmip2_cfc.o ocmip2_biotic.o ocean_bgc_restore.o ocmip2_he.o ocean_po4_pre.o ocean_ibgc.o ocean_generic_tracer.o ocean_frazil.o ocean_tempsalt.o ocean_passive.o transport_matrix.o time_manager.o coupler_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_tpm.F90
ocean_tpm_util.o: ./ocean_tpm_util.F90 field_manager.o fms.o mpp.o fm_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_tpm_util.F90
ocean_tracer.o: ./ocean_tracer.F90 constants.o diag_manager.o field_manager.o fms.o fms_io.o fm_util.o mpp_domains.o mpp_io.o mpp.o platform.o time_manager.o transport_matrix.o ocean_blob.o ocean_convect.o ocean_domains.o ocean_density.o ocean_frazil.o ocean_bih_tracer.o ocean_generic_tracer.o ocean_lap_tracer.o ocean_obc.o ocean_parameters.o ocean_passive.o ocean_shortwave.o ocean_tempsalt.o ocean_tpm.o ocean_tpm_util.o ocean_tracer_advect.o ocean_tracer_diag.o ocean_tracer_util.o ocean_thickness.o ocean_types.o ocean_util.o ocean_vert_mix.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_tracer.F90
ocean_tracer_advect.o: ./ocean_tracer_advect.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_domains.o ocean_obc.o ocean_parameters.o ocean_topog.o ocean_tracer_util.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_tracer_advect.F90
ocean_tracer_diag.o: ./ocean_tracer_diag.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o time_manager.o ocean_density.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_tracer_util.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_tracer_diag.F90
ocean_tracer_util.o: ./ocean_tracer_util.F90 constants.o fms.o mpp.o platform.o diag_manager.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_tracer_util.F90
ocean_types.o: ./ocean_types.F90 coupler_types.o field_manager.o fms.o mpp_domains.o mpp.o time_manager.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_types.F90
ocean_util.o: ./ocean_util.F90 platform.o constants.o diag_manager.o mpp_domains.o mpp.o time_manager.o ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_util.F90
ocean_velocity.o: ./ocean_velocity.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o time_interp_external.o ocean_bih_friction.o ocean_coriolis.o ocean_domains.o ocean_form_drag.o ocean_lap_friction.o ocean_momentum_source.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_pressure.o ocean_types.o ocean_util.o ocean_velocity_advect.o ocean_velocity_diag.o ocean_vert_mix.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_velocity.F90
ocean_velocity_advect.o: ./ocean_velocity_advect.F90 diag_manager.o fms.o mpp.o ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_velocity_advect.F90
ocean_velocity_diag.o: ./ocean_velocity_diag.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o time_manager.o ocean_bih_friction.o ocean_coriolis.o ocean_domains.o ocean_form_drag.o ocean_lap_friction.o ocean_momentum_source.o ocean_operators.o ocean_parameters.o ocean_pressure.o ocean_types.o ocean_util.o ocean_velocity_advect.o ocean_vert_mix.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_velocity_diag.F90
ocean_vert_chen.o: ./ocean_vert_chen.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_shortwave_csiro.o ocean_types.o ocean_util.o ocean_vert_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_vert_chen.F90
ocean_vert_const.o: ./ocean_vert_const.F90 constants.o diag_manager.o fms.o mpp_io.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_vert_const.F90
ocean_vert_gotm.o: ./ocean_vert_gotm.F90 constants.o diag_manager.o fms.o fms_io.o mpp_domains.o mpp.o tridiagonal_gotm.o turbulence.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o ocean_obc.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_vert_gotm.F90
ocean_vert_kpp_mom4p0.o: ./ocean_vert_kpp_mom4p0.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_vert_kpp_mom4p0.F90
ocean_vert_kpp_mom4p1.o: ./ocean_vert_kpp_mom4p1.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_vert_kpp_mom4p1.F90
ocean_vert_kpp_test.o: ./ocean_vert_kpp_test.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_vert_util.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_vert_kpp_test.F90
ocean_vert_mix.o: ./ocean_vert_mix.F90 constants.o diag_manager.o field_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_parameters.o ocean_tracer_util.o ocean_types.o ocean_util.o ocean_vert_const.o ocean_vert_chen.o ocean_vert_gotm.o ocean_vert_kpp_test.o ocean_vert_kpp_mom4p0.o ocean_vert_kpp_mom4p1.o ocean_vert_pp.o ocean_vert_tidal.o ocean_vert_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_vert_mix.F90
ocean_vert_pp.o: ./ocean_vert_pp.F90 constants.o diag_manager.o fms.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_vert_util.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_vert_pp.F90
ocean_vert_tidal.o: ./ocean_vert_tidal.F90 constants.o diag_manager.o fms.o mpp_domains.o mpp.o ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_vert_tidal.F90
ocean_vert_util.o: ./ocean_vert_util.F90 constants.o diag_manager.o fms.o mpp.o ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_vert_util.F90
ocean_workspace.o: ./ocean_workspace.F90 ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocean_workspace.F90
ocmip2_abiotic.o: ./ocmip2_abiotic.F90 atmos_ocean_fluxes.o diag_manager.o field_manager.o time_manager.o fms.o fms_io.o mpp.o time_interp_external.o mpp_domains.o constants.o ocean_tpm_util.o fm_util.o coupler_types.o ocean_types.o ocmip2_co2calc.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocmip2_abiotic.F90
ocmip2_biotic.o: ./ocmip2_biotic.F90 time_manager.o diag_manager.o field_manager.o fms.o fms_io.o mpp.o time_interp_external.o mpp_domains.o constants.o ocean_tpm_util.o fm_util.o coupler_types.o ocean_types.o ocmip2_co2calc.o ocean_util.o atmos_ocean_fluxes.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocmip2_biotic.F90
ocmip2_cfc.o: ./ocmip2_cfc.F90 atmos_ocean_fluxes.o time_manager.o field_manager.o fms.o ocean_tpm_util.o constants.o fm_util.o mpp.o coupler_types.o ocean_types.o diag_manager.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocmip2_cfc.F90
ocmip2_co2calc.o: ./ocmip2_co2calc.F90
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocmip2_co2calc.F90
ocmip2_he.o: ./ocmip2_he.F90 atmos_ocean_fluxes.o time_manager.o time_interp_external.o diag_manager.o field_manager.o fms.o fms_io.o ocean_tpm_util.o fm_util.o mpp.o coupler_types.o ocean_types.o mpp_domains.o ocean_util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./ocmip2_he.F90
platform.o: ./platform.F90
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./platform.F90
potentialml.o: ./potentialml.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./potentialml.F90
production.o: ./production.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./production.F90
q2over2eq.o: ./q2over2eq.F90 turbulence.o util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./q2over2eq.F90
read_mosaic.o: ./read_mosaic.c
	$(CC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	./read_mosaic.c
threadloc.o: ./threadloc.c
	$(CC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	./threadloc.c
time_interp.o: ./time_interp.F90 time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./time_interp.F90
time_interp_external.o: ./time_interp_external.F90 mpp.o mpp_io.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils.o fms.o platform.o horiz_interp.o constants.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./time_interp_external.F90
time_manager.o: ./time_manager.F90 constants.o fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./time_manager.F90
tkealgebraic.o: ./tkealgebraic.F90 turbulence.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./tkealgebraic.F90
tkeeq.o: ./tkeeq.F90 turbulence.o util.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./tkeeq.F90
tracer_manager.o: ./tracer_manager.F90 mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./tracer_manager.F90
transport_matrix.o: ./transport_matrix.F90 field_manager.o mpp.o diag_manager.o ocean_tpm_util.o ocean_types.o
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./transport_matrix.F90
tridiagonal_gotm.o: ./tridiagonal_gotm.F90
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./tridiagonal_gotm.F90
turbulence.o: ./turbulence.F90
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./turbulence.F90
util.o: ./util.F90
	$(FC) $(CPPDEFS) $(TOPAZ_CPPFLAGS) $(FPPFLAGS) $(TOPAZ_FFLAGS) $(OTHERFLAGS) -c	./util.F90
SRC = ./ocean_residency_range.F90 ./ocean_blob_dynamic_free.F90 ./ispralength.F90 ./fms_io.F90 ./time_interp_external.F90 ./ocean_vert_gotm.F90 ./ocean_residency_integrand.F90 ./coupler_types.F90 ./ocmip2_biotic.F90 ./mpp_memutils.F90 ./ocean_blob_diag.F90 ./diag_util.F90 ./ocean_velocity_diag.F90 ./FMS_coupler_util.F90 ./tkeeq.F90 ./generic_TOPAZ_model.F90 ./ocean_shortwave_gfdl.F90 ./time_interp.F90 ./genericeq.F90 ./mosaic.F90 ./cmue_rf.F90 ./memutils.F90 ./ocean_lap_tracer.F90 ./horiz_interp_bicubic.F90 ./ocean_vert_kpp_mom4p1.F90 ./kbalgebraic.F90 ./util.F90 ./diag_axis.F90 ./diag_manager.F90 ./ocean_lapcst_friction.F90 ./ocean_residency.F90 ./ocean_form_drag.F90 ./ocean_vert_chen.F90 ./create_xgrid.c ./ocean_tempsalt.F90 ./horiz_interp_conserve.F90 ./kbeq.F90 ./ocean_age_tracer.F90 ./generic_ERGOM.F90 ./mpp_parameter.F90 ./axis_utils.F90 ./atmos_ocean_fluxes.F90 ./platform.F90 ./netcdf_var.F90 ./mpp_domains.F90 ./cmue_c.F90 ./ocean_thickness.F90 ./ocean_bih_tracer.F90 ./ocean_vert_kpp_test.F90 ./ocean_velocity_advect.F90 ./diag_grid.F90 ./ocean_tracer_diag.F90 ./ocean_vert_const.F90 ./ocean_residency_meta.F90 ./ocean_density.F90 ./alpha_mnb.F90 ./ocean_workspace.F90 ./ocean_operators.F90 ./FMS_ocmip2_co2calc.F90 ./potentialml.F90 ./mpp_io.F90 ./ocean_coriolis.F90 ./ocean_bihcgrid_friction.F90 ./fm_util.F90 ./ocean_momentum_source.F90 ./diag_data.F90 ./read_mosaic.c ./generic_tracer_utils.F90 ./get_cal_time.F90 ./generic_CFC.F90 ./ocean_pert_co2.F90 ./constants.F90 ./ocean_pressure.F90 ./ocean_lapgen_friction.F90 ./generic_TOPAZ.F90 ./ocean_vert_util.F90 ./internal_wave.F90 ./production.F90 ./cmue_d.F90 ./tkealgebraic.F90 ./ocean_tpm_util.F90 ./ocean_lapcgrid_friction.F90 ./mpp_data.F90 ./ocean_passive.F90 ./compute_rist.F90 ./data_override.F90 ./ocean_obc.F90 ./tracer_manager.F90 ./ocean_bgc_restore.F90 ./ocean_vert_mix.F90 ./mosaic_util.c ./ocean_domains.F90 ./ocean_shortwave.F90 ./ocean_velocity.F90 ./ocean_grids.F90 ./ocean_topog.F90 ./ocmip2_abiotic.F90 ./lengthscaleeq.F90 ./diff_face.F90 ./ocean_shortwave_csiro.F90 ./threadloc.c ./ocmip2_he.F90 ./ocmip2_cfc.F90 ./ocean_generic_tracer.F90 ./ocean_lap_friction.F90 ./generic_BLING.F90 ./ocean_parameters.F90 ./ocean_util.F90 ./ocean_vert_pp.F90 ./ocean_shortwave_jerlov.F90 ./ocean_vert_kpp_mom4p0.F90 ./ocean_blob_dynamic_bottom.F90 ./tridiagonal_gotm.F90 ./fms.F90 ./mpp_pset.F90 ./compute_cpsi3.F90 ./ocean_blob_util.F90 ./diag_output.F90 ./ocean_tracer.F90 ./ocean_residency_ml.F90 ./ocean_obc_barotrop.F90 ./dissipationeq.F90 ./ocean_convect.F90 ./diag_table.F90 ./time_manager.F90 ./cmue_sg.F90 ./ocean_blob_static_free.F90 ./ocean_bihcst_friction.F90 ./horiz_interp_spherical.F90 ./netcdf_io.F90 ./ocean_tracer_util.F90 ./ocean_barotropic.F90 ./ocean_bih_friction.F90 ./horiz_interp.F90 ./ocean_bihgen_friction.F90 ./ocean_frazil.F90 ./ocean_po4_pre.F90 ./ocean_tpm.F90 ./algebraiclength.F90 ./ocean_blob.F90 ./ocean_blob_static_bottom.F90 ./horiz_interp_bilinear.F90 ./generic_tracer.F90 ./cmue_ma.F90 ./epsbalgebraic.F90 ./ocean_tracer_advect.F90 ./ocmip2_co2calc.F90 ./ocean_vert_tidal.F90 ./field_manager.F90 ./ocean_types.F90 ./turbulence.F90 ./q2over2eq.F90 ./transport_matrix.F90 ./ocean_ibgc.F90 ./mpp.F90 ./fk_craig.F90 ./horiz_interp_type.F90
OBJ = ocean_residency_range.o ocean_blob_dynamic_free.o ispralength.o fms_io.o time_interp_external.o ocean_vert_gotm.o ocean_residency_integrand.o coupler_types.o ocmip2_biotic.o mpp_memutils.o ocean_blob_diag.o diag_util.o ocean_velocity_diag.o FMS_coupler_util.o tkeeq.o generic_TOPAZ_model.o ocean_shortwave_gfdl.o time_interp.o genericeq.o mosaic.o cmue_rf.o memutils.o ocean_lap_tracer.o horiz_interp_bicubic.o ocean_vert_kpp_mom4p1.o kbalgebraic.o util.o diag_axis.o diag_manager.o ocean_lapcst_friction.o ocean_residency.o ocean_form_drag.o ocean_vert_chen.o create_xgrid.o ocean_tempsalt.o horiz_interp_conserve.o kbeq.o ocean_age_tracer.o generic_ERGOM.o mpp_parameter.o axis_utils.o atmos_ocean_fluxes.o platform.o netcdf_var.o mpp_domains.o cmue_c.o ocean_thickness.o ocean_bih_tracer.o ocean_vert_kpp_test.o ocean_velocity_advect.o diag_grid.o ocean_tracer_diag.o ocean_vert_const.o ocean_residency_meta.o ocean_density.o alpha_mnb.o ocean_workspace.o ocean_operators.o FMS_ocmip2_co2calc.o potentialml.o mpp_io.o ocean_coriolis.o ocean_bihcgrid_friction.o fm_util.o ocean_momentum_source.o diag_data.o read_mosaic.o generic_tracer_utils.o get_cal_time.o generic_CFC.o ocean_pert_co2.o constants.o ocean_pressure.o ocean_lapgen_friction.o generic_TOPAZ.o ocean_vert_util.o internal_wave.o production.o cmue_d.o tkealgebraic.o ocean_tpm_util.o ocean_lapcgrid_friction.o mpp_data.o ocean_passive.o compute_rist.o data_override.o ocean_obc.o tracer_manager.o ocean_bgc_restore.o ocean_vert_mix.o mosaic_util.o ocean_domains.o ocean_shortwave.o ocean_velocity.o ocean_grids.o ocean_topog.o ocmip2_abiotic.o lengthscaleeq.o diff_face.o ocean_shortwave_csiro.o threadloc.o ocmip2_he.o ocmip2_cfc.o ocean_generic_tracer.o ocean_lap_friction.o generic_BLING.o ocean_parameters.o ocean_util.o ocean_vert_pp.o ocean_shortwave_jerlov.o ocean_vert_kpp_mom4p0.o ocean_blob_dynamic_bottom.o tridiagonal_gotm.o fms.o mpp_pset.o compute_cpsi3.o ocean_blob_util.o diag_output.o ocean_tracer.o ocean_residency_ml.o ocean_obc_barotrop.o dissipationeq.o ocean_convect.o diag_table.o time_manager.o cmue_sg.o ocean_blob_static_free.o ocean_bihcst_friction.o horiz_interp_spherical.o netcdf_io.o ocean_tracer_util.o ocean_barotropic.o ocean_bih_friction.o horiz_interp.o ocean_bihgen_friction.o ocean_frazil.o ocean_po4_pre.o ocean_tpm.o algebraiclength.o ocean_blob.o ocean_blob_static_bottom.o horiz_interp_bilinear.o generic_tracer.o cmue_ma.o epsbalgebraic.o ocean_tracer_advect.o ocmip2_co2calc.o ocean_vert_tidal.o field_manager.o ocean_types.o turbulence.o q2over2eq.o transport_matrix.o ocean_ibgc.o mpp.o fk_craig.o horiz_interp_type.o
realclean: neat
	-rm -f .libtopaz.a.cppdefs $(OBJ) libtopaz.a *.mod
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
libtopaz.a: $(OBJ) $(LIBDIR)/libbio_prod.a $(LIBDIR)/libturbulence_prod.a $(LIBDIR)/libutil_prod.a
	$(AR) $(ARFLAGS) libtopaz.a $(OBJ) $(LIBDIR)/libbio_prod.a $(LIBDIR)/libturbulence_prod.a $(LIBDIR)/libutil_prod.a
	cp *.mod $(MODDIR)
	cp libtopaz.a $(LIBDIR)
