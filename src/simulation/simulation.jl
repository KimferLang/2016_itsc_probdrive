export  SimParams,
		DEFAULT_SIM_PARAMS,

        simulate!,
        simulate_but_terminate_if_collision!,
        propagate!,

        get_input_acceleration,
		get_input_turnrate,

        calc_sequential_moving_average,
		calc_weighted_moving_average

immutable SimParams
	sec_per_frame :: Float64
	n_euler_steps :: Int

	function SimParams(
		sec_per_frame :: Float64 = DEFAULT_SEC_PER_FRAME,
		n_euler_steps :: Int     = 10,
		)

		@assert(sec_per_frame > 0.0)
		@assert(n_euler_steps > 0)

		new(sec_per_frame, n_euler_steps)
	end
end
const DEFAULT_SIM_PARAMS = SimParams()

function simulate!(
    runlog          :: RunLog,
    sn              :: StreetNetwork,
    behavior        :: AbstractVehicleBehavior,
    id              :: UInt,
    frame_start     :: Int,
    frame_end       :: Int;
    pdset_frames_per_sim_frame::Int = N_FRAMES_PER_SIM_FRAME,
    n_euler_steps   :: Int = 2
    )

    targets = get_targets(behavior)

    for frame in frame_start : pdset_frames_per_sim_frame : frame_end-1

        colset = id2colset(runlog, id, frame)
        @assert(colset != COLSET_NULL)

        action_lat, action_lon = select_action(behavior, runlog, sn, colset, frame)
        propagate!(runlog, sn, frame, colset, targets, action_lat, action_lon,
                   pdset_frames_per_sim_frame, n_euler_steps)
    end

    runlog
end

get_input_acceleration(::Features.Feature_FutureAcceleration, action_lon::Float64, acc_prev::Float64) = action_lon
get_input_acceleration(::Features.Feature_Future_Delta_Accel, action_lon::Float64, acc_prev::Float64) = acc_prev + action_lon

get_input_turnrate(::Features.Feature_FutureTurnrate, action_lat::Float64, ϕ::Float64, v::Float64, a::Float64, δt::Float64) = action_lat
function get_input_turnrate(::Features.Feature_FutureVelFt, action_lat::Float64, ϕ::Float64, v::Float64, a::Float64, δt::Float64)
    v_lat = action_lat
    (asin(v_lat/(v + a*δt)) - ϕ)/δt
end
function get_input_turnrate(::Features.Feature_FutureDesiredAngle, action_lat::Float64, ϕ::Float64, v::Float64, a::Float64, δt::Float64)
    phi_des = action_lat
    (phi_des - ϕ)*Features.KP_DESIRED_ANGLE
end

function _propagate_one_runlog_frame!(
    runlog        :: RunLog,
    sn            :: StreetNetwork,
    frame         :: Int,
    colset        :: UInt,
    targets       :: ModelTargets,
    action_lat    :: Float64,
    action_lon    :: Float64,
    n_euler_steps :: Int,
    )

    frame_fut = frame + 1
    Δt = RunLogs.get_elapsed_time(runlog, frame, frame_fut)
    δt = Δt / n_euler_steps

    inertial = get(runlog, colset, frame, :inertial)::VecSE2
    x = inertial.x
    y = inertial.y
    θ = inertial.θ

    s = d = 0.0
    ϕ = ϕₒ = (get(runlog, colset, frame, :frenet)::VecSE2).θ

    rates = get(runlog, colset, frame, :ratesB)::VecSE2
    v = sqrt(rates.x*rates.x + rates.y*rates.y)

    acc_prev = (hypot(get(runlog, id2colset(runlog, colset2id(runlog, colset, frame), frame-1), frame-1, :ratesB)) - v)/Δt
    for i = 1 : n_euler_steps

        a = get_input_acceleration(targets.lon, action_lon, acc_prev)
        ω = get_input_turnrate(targets.lat, action_lat, ϕ, v, a, δt)

        v += a*δt
        θ += ω*δt
        x += v*cos(θ)*δt
        y += v*sin(θ)*δt

        proj = project_point_to_streetmap(x, y, sn)
        @assert(proj.successful)
        s, d, ϕ = pt_to_frenet_xyy(proj.footpoint, x, y, θ)
        acc_prev = a
    end

    proj = project_point_to_streetmap(x, y, sn)
    @assert(proj.successful)
    s, d, ϕ = pt_to_frenet_xyy(proj.footpoint, x, y, θ)

    colset_fut = get(runlog, colset, frame, :next_colset)::UInt
    if colset_fut == COLSET_NULL
        # automatically insert car into future frame if necessary
        id = colset2id(runlog, colset, frame)
        colset_fut = get_first_vacant_colset!(runlog, id, frame_fut)
    end

    RunLogs.set!(runlog, colset_fut, frame_fut,
         get(runlog, colset, frame, :id)::UInt,
         VecSE2(x, y, θ), # inertial
         VecSE2(s, d, ϕ), # frenet
         VecSE2(v, 0.0, (ϕ - ϕₒ)/Δt), # ratesB (assume zero sideslip)
         proj.extind,
         proj.footpoint,
         proj.lane.id,
         COLSET_NULL,
         COLSET_NULL,
         get(runlog, colset, frame, :behavior)::UInt16
    )

    RunLogs.set!(runlog, colset_fut, frame_fut, :colset_front,
        RunLogs.calc_front_vehicle_colset(runlog, sn, colset_fut, frame_fut))
    RunLogs.set!(runlog, colset_fut, frame_fut, :colset_rear,
        RunLogs.calc_rear_vehicle_colset(runlog, sn, colset_fut, frame_fut))

    runlog
end
function propagate!(
    runlog        :: RunLog,
    sn            :: StreetNetwork,
    frame         :: Int,
    colset        :: UInt,
    targets       :: ModelTargets,
    action_lat    :: Float64,
    action_lon    :: Float64,
    pdset_frames_per_sim_frame :: Int,
    n_euler_steps :: Int,
    )

    for jump in 0 : pdset_frames_per_sim_frame-1
        frame_fut = frame + jump
        _propagate_one_runlog_frame!(runlog, sn, frame_fut, colset, targets, action_lat, action_lon, n_euler_steps)
    end

    runlog
end

# TODO(tim): move to elsewhere?
function calc_sequential_moving_average(
	vec         :: AbstractArray{Float64}, # vector of values to smooth on
	index_start :: Int,                    # the present index; value must be already populated
	history     :: Int                     # the number of values to smooth over, (≥ 1)
	)

	# Sequential Moving Average: the average of the past n results

	@assert(history ≥ 1)

	clamped_history = min(history, index_start)
	index_low = index_start - clamped_history + 1

	retval = 0.0
	for i = index_low : index_start
		retval += vec[i]
	end
	retval / clamped_history
end
function calc_weighted_moving_average(
	vec         :: AbstractArray{Float64}, # vector of values to smooth on
	index_start :: Int,                    # the present index; value must be already populated
	history     :: Int                     # the number of values to smooth over, (≥ 1)
	)

	# Weighted Moving Average: the average of the past n results weighted linearly
	# ex: (3×f₁ + 2×f₂ + 1×f₃) / (3 + 2 + 1)

	@assert(history ≥ 1)

	clamped_history = min(history, index_start)
	index_low = index_start - clamped_history + 1

	retval = 0.0
	for i = index_low : index_start
		retval += vec[i] * (i - index_low + 1)
	end
	retval / (0.5clamped_history*(clamped_history+1))
end

function _reverse_smoothing_sequential_moving_average(
	vec::AbstractArray{Float64}, # vector of values originally smoothed on;
	                             # with the most recent value having been overwritten with the smoothed value
	index_start::Int, # the present index; value must be already populated
	history::Int # the number of values to smooth over, (≥ 1)
	)

	# If the SMA is (f₁ + f₂ + f₃ + ...) / n
	# the reverse value is f₁ = n⋅SMA - f₂ - f₃ - ...

	@assert(history ≥ 1)

	clamped_history = min(history, index_start)
	index_low = index_start - clamped_history + 1

	smoothed_result = vec[index_start]

	retval = clamped_history * smoothed_result
	for i = index_low : index_start-1
		retval -= vec[i]
	end
	retval
end
function _reverse_smoothing_weighted_moving_average(
	vec::AbstractArray{Float64}, # vector of values originally smoothed on;
	                             # with the most recent value having been overwritten with the smoothed value
	index_start::Int, # the present index; value must be already populated
	history::Int # the number of values to smooth over, (≥ 1)
	)

	# If the WMA is (3×f₁ + 2×f₂ + 1×f₃) / (3 + 2 + 1)
	# the reverse value is f₁ = [WMA * (3 + 2 + 1) - 2×f₂ - 1×f₃] / 3

	@assert(history ≥ 1)

	clamped_history = min(history, index_start)
	index_low = index_start - clamped_history + 1

	smoothed_result = vec[index_start]

	retval = (0.5clamped_history*(clamped_history+1)) * smoothed_result
	for i = index_low : index_start-1
		retval -= vec[i] * (i - index_low + 1)
	end
	retval / clamped_history
end

