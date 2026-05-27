function effectiveSchedule(::AbstractInference, schedule)
    return schedule === nothing ? :sequential : schedule
end

function validateGBPSchedule(schedule)
    if !(schedule in (:sequential, :flooding, :residual))
        error("Schedule must be :sequential, :flooding, or :residual.")
    end

    return nothing
end