module modRotation
  export rotationFactory,
         rotator_Yes,
         rotator_No

  abstract type Rotator end

  struct rotator_Yes <: Rotator end

  struct rotator_No  <: Rotator end

  function rotationFactory(type::String)
    if type == "Yes"
      return rotator_Yes()
    else
      return rotator_No()
    end
  end

end