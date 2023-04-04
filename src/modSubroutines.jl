module modSubroutines
  export get_tags

  using Gridap.Geometry

  function get_tags(matFlag, labels, dimens)
    tags = copy(get_face_tag(labels,dimens))
    if length(matFlag) != 0
      for (iflag, flag) in enumerate(matFlag)
          num_tag = get_tag_from_name(labels,flag)
          replace!(tags, num_tag => iflag)
      end
    else
      tags .= 1
    end
    return tags
  end

end